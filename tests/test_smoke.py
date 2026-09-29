"""Offline runtime checks, not validation of the scientific model.

Run from the repository root: python -m unittest discover -s tests -v -b
The expected failure records the unresolved matrix-indexing issue described
in KNOWN_ISSUES.md; it must be revisited when that issue is corrected.
"""

import io
import json
import math
import os
from pathlib import Path
import re
import tempfile
import unittest


class WebSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temporary = tempfile.TemporaryDirectory(prefix="helix-tests-")
        cls.addClassCleanup(cls.temporary.cleanup)
        cls.previous_env = {
            key: os.environ.get(key)
            for key in ("HELIX_CACHE_DIR", "HELIX_VAR_DIR", "HELIX_USAGE_LOGGING")
        }
        os.environ["HELIX_CACHE_DIR"] = str(Path(cls.temporary.name) / "cache")
        os.environ["HELIX_VAR_DIR"] = str(Path(cls.temporary.name) / "var")
        os.environ["HELIX_USAGE_LOGGING"] = "0"
        cls.addClassCleanup(cls.restore_environment)
        # Configuration is read at import time, after setting temporary paths.
        from fastapi.testclient import TestClient
        from app.main import app
        from app import engine

        cls.engine = engine
        cls.client = TestClient(app)
        cls.addClassCleanup(cls.client.close)
        cls.common = {"residues": engine.DEFAULT_ALPHABET_16, "n_decoys": 100, "seed": 42}
        cls.design = {**cls.common, "length": 4, "steps": 1000, "restarts": 1}
        cls.env_a = {"preset": "interfacial_neg"}
        cls.env_b = {"preset": "homogeneous_polar"}

    @classmethod
    def restore_environment(cls):
        for key, value in cls.previous_env.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value

    def assert_finite(self, value):
        if isinstance(value, dict):
            for child in value.values():
                self.assert_finite(child)
        elif isinstance(value, list):
            for child in value:
                self.assert_finite(child)
        elif isinstance(value, (int, float)):
            self.assertTrue(math.isfinite(value), value)

    def post(self, route, payload):
        response = self.client.post(route, json=payload)
        self.assertEqual(response.status_code, 200, response.text)
        result = response.json()
        self.assert_finite(result)
        return result

    def stream(self, route, payload):
        with self.client.stream("POST", route, json=payload) as response:
            self.assertEqual(response.status_code, 200)
            self.assertIn("application/x-ndjson", response.headers["content-type"])
            events = [json.loads(line) for line in response.iter_lines() if line]
        self.assertTrue(events)
        self.assertFalse([event for event in events if event["type"] == "error"], events)
        self.assertEqual(events[-1]["type"], "final")
        self.assertEqual(sum(event["type"] == "final" for event in events), 1)
        self.assert_finite(events)
        return events[-1]["payload"]

    def assert_score(self, result):
        self.assertEqual(len(result["sequence"]), 4)
        self.assertEqual(result["length"], 4)
        self.assertTrue(set(result["sequence"]) <= set(self.common["residues"]))
        self.assertAlmostEqual(result["energy"], result["breakdown"]["total"])
        self.assertGreater(result["random_reference"]["std"], 0)
        self.assertLessEqual(0, result["percentile_lower_is_better"])
        self.assertLessEqual(result["percentile_lower_is_better"], 100)
        self.assertEqual(len(result["wheel"]), 4)

    def test_page_assets_and_health(self):
        response = self.client.get("/")
        self.assertEqual(response.status_code, 200)
        assets = set(re.findall(r'(?:src|href)="(/static/[^"?#]+)', response.text))
        self.assertTrue(assets)
        for asset in assets:
            with self.subTest(asset=asset):
                response = self.client.get(asset)
                self.assertEqual(response.status_code, 200)
                self.assertTrue(response.content)
        health = self.client.get("/api/health")
        self.assertEqual(health.status_code, 200)
        self.assertEqual(health.json()["status"], "ok")
        self.assertEqual(self.client.get("/openapi.json").status_code, 200)

    def test_score_and_reproducibility(self):
        payload = {**self.common, "sequence": "ALKD", "n_random": 100}
        first = self.post("/api/score", payload)
        self.assert_score(first)
        self.assertEqual(first, self.post("/api/score", payload))

    def test_single_design(self):
        self.assert_score(self.post("/api/design", self.design))

    def test_comparison_and_cross_design(self):
        environments = {"environment_a": self.env_a, "environment_b": self.env_b}
        requests = {
            "/api/compare": {**self.common, **environments, "sequence": "ALKD"},
            "/api/cross-design": {**self.design, **environments},
        }
        for route, payload in requests.items():
            with self.subTest(route=route):
                result = self.post(route, payload)
                self.assert_score(result["environment_a"])
                self.assert_score(result["environment_b"])
                self.assertEqual(result["environment_a"]["sequence"], result["sequence"])
                self.assertEqual(result["environment_b"]["sequence"], result["sequence"])

    def test_family_and_stream(self):
        payload = {**self.design, "family_size": 2, "oversample_factor": 1}
        result = self.post("/api/design-family", payload)
        self.assertEqual(result["family_size_returned"], 2)
        self.assertEqual(len({member["sequence"] for member in result["members"]}), 2)
        for member in result["members"]:
            self.assert_score(member)
        self.assertEqual(result, self.stream("/api/design-family-stream", payload))

    def test_specificity_and_stream(self):
        payload = {
            **self.design,
            "target_environment": self.env_a,
            "off_target_environments": [self.env_b],
            "num_sequences": 1,
        }
        result = self.post("/api/specificity-design", payload)
        self.assertEqual(result["specificity_collection"]["returned"], 1)
        self.assertEqual(len(result["candidates"]), 1)
        self.assert_score(result["target_environment"])
        self.assertEqual(len(result["off_target_environments"]), 1)
        self.assert_score(result["off_target_environments"][0])
        self.assertEqual(result, self.stream("/api/specificity-design-stream", payload))

    def test_penetration_scan(self):
        result = self.post("/api/optimal-penetration", {
            **self.common, "sequence": "ALKD", "environment": self.env_a, "n_random": 100,
        })
        self.assertEqual([point["penetration_percent"] for point in result["profile"]], list(range(0, 101, 5)))
        self.assert_score(result["best_score"])
        best = min(result["profile"], key=lambda point: point["z_score"])
        self.assertEqual(result["best_penetration_percent"], best["penetration_percent"])

    def test_pdb_export(self):
        from Bio.PDB import PDBParser

        result = self.post("/api/pdb", {"sequence": "ALKD"})
        structure = PDBParser(QUIET=True).get_structure("test", io.StringIO(result["pdb"]))
        self.assertEqual([residue.resname for residue in structure.get_residues()], ["ALA", "LEU", "LYS", "ASP"])
        for atom in structure.get_atoms():
            self.assertTrue(all(math.isfinite(float(x)) for x in atom.coord))

    def test_request_limits(self):
        invalid = [
            ("/api/score", {"sequence": ""}),
            ("/api/design", {"length": 41}),
            ("/api/design", {"residues": ["A"]}),
            ("/api/optimal-penetration", {"sequence": "ALKD", "environment": self.env_b}),
        ]
        for route, payload in invalid:
            with self.subTest(route=route, payload=payload):
                self.assertEqual(self.client.post(route, json=payload).status_code, 422)

    def test_usage_logging_disabled(self):
        self.assertEqual(self.client.post("/api/visit").status_code, 200)
        self.assertEqual(self.client.get("/api/usage-metrics").status_code, 200)
        self.assertFalse((Path(self.temporary.name) / "var" / "usage_events.jsonl").exists())

    @unittest.expectedFailure
    def test_known_neighbor_score_alphabet_order_invariance(self):
        """Raw chemical scores should not change when only alphabet order changes."""
        alphabet = self.engine.DEFAULT_ALPHABET_16
        reordered = ["P"] + [aa for aa in alphabet if aa != "P"]
        env = self.engine.resolve_environment(self.env_a)
        values = []
        for residues in (alphabet, reordered):
            builder = self.engine.build_builder(4, residues, env, n_decoys=100)
            values.append(builder._raw_helix_neigh(self.engine.seq_to_codes("PPPP", residues)))
        self.assertAlmostEqual(values[0], values[1])


if __name__ == "__main__":
    unittest.main()
