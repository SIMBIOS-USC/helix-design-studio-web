"""Search and displayed scores must describe the same calibrated objective."""

import tempfile
import unittest


class SearchConsistencyTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory(prefix="helix-search-test-")
        self.addCleanup(self.temporary.cleanup)
        from app import engine
        self.engine = engine
        # The engine may already be imported by other tests; override the one
        # default path captured when its builder factory was defined.
        defaults = engine.build_builder.__defaults__
        engine.build_builder.__defaults__ = (*defaults[:-1], self.temporary.name)
        self.addCleanup(setattr, engine.build_builder, "__defaults__", defaults)
        self.settings = {
            "length": 6, "residues": list("AEKLP"), "n_decoys": 100,
            "steps": 1000, "restarts": 2, "seed": 71,
        }
        self.target = {"preset": "interfacial_neg"}
        self.off_target = {"preset": "homogeneous_polar"}

    def test_cross_objective_matches_displayed_components(self):
        result = self.engine.cross_environment_design(
            **self.settings, environment_a=self.target,
            environment_b=self.off_target, lambda_gap=0.75,
        )
        a = result["environment_a"]["energy"]
        b = result["environment_b"]["energy"]
        self.assertAlmostEqual(result["cross_design"]["energy_a"], a, places=12)
        self.assertAlmostEqual(result["cross_design"]["energy_b"], b, places=12)
        self.assertAlmostEqual(result["cross_design"]["objective"], a+b+0.75*abs(a-b), places=12)
        self.assertEqual(result["environment_b"]["calibration"]["seed"], 71)

    def test_specificity_objective_matches_every_candidate(self):
        result = self.engine.design_specificity(
            **self.settings, target_environment=self.target,
            off_target_environments=[self.off_target],
            lambda_balance=0.4, num_sequences=2,
        )
        for candidate in result["candidates"]:
            with self.subTest(sequence=candidate["sequence"]):
                target = candidate["target_environment"]
                competitors = candidate["off_target_environments"]
                worst = min(score["energy"] for score in competitors)
                details = candidate["specificity_design"]
                self.assertAlmostEqual(details["best_target_energy"], target["energy"], places=12)
                self.assertAlmostEqual(details["best_off_target_energy"], worst, places=12)
                self.assertAlmostEqual(details["objective"], target["energy"]-0.6*worst, places=12)
                self.assertAlmostEqual(details["energy_margin"], worst-target["energy"], places=12)
                self.assertEqual(target["orientation_mode"], "fixed")
                self.assertEqual(target["calibration"]["n_decoys"], 100)
                self.assertEqual(target["calibration"]["seed"], 71)


if __name__ == "__main__":
    unittest.main()
