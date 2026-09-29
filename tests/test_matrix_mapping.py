"""Regression checks for parameter indexing, not scientific validation.

M3/M4 expectations are conditional on the provisional header reconstruction
documented in PARAMETER_PROVENANCE.md. No test establishes its provenance.
"""

from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "Code"))

from core.hamiltonian_builder import HamiltonianBuilder
from data_loaders import energy_matrix_loader as loader


FULL_ALPHABET = list("ACDEFGHIKLMNPQRSTVWY")
DEFAULT_ALPHABET = list("APDEQLMFKSTWYVHG")


def make_builder(alphabet, length=5):
    # Indexing tests do not need Monte Carlo calibration or persistent caches.
    with patch.object(HamiltonianBuilder, "_load_or_calibrate_z_scores", return_value={}):
        return HamiltonianBuilder(length, list(alphabet), 5, 5 * length)


class MatrixMappingTests(unittest.TestCase):
    def test_all_matrix_entries_follow_residue_identity(self):
        loaders = {
            "mj_matrix": loader._load_energy_matrix_file,
            "M1": loader._load_first_neighbors_matrix_file,
            "M3": loader._load_third_neighbors_matrix_file,
            "M4": loader._load_fourth_neighbors_matrix_file,
        }
        alphabets = (
            FULL_ALPHABET,
            FULL_ALPHABET[::-1],
            FULL_ALPHABET[7:] + FULL_ALPHABET[:7],
            DEFAULT_ALPHABET,
            list("PCYAGL"),
        )
        for alphabet in alphabets:
            builder = make_builder(alphabet)
            for name, load in loaders.items():
                matrix, symbols = load()
                lookup = {
                    (aa, bb): matrix[i, j]
                    for i, aa in enumerate(symbols)
                    for j, bb in enumerate(symbols)
                }
                actual = getattr(builder, name)
                with self.subTest(alphabet=alphabet, matrix=name):
                    self.assertEqual(actual.shape, (len(alphabet), len(alphabet)))
                    for i, aa in enumerate(alphabet):
                        for j, bb in enumerate(alphabet):
                            self.assertEqual(actual[i, j], lookup[aa, bb])

    def test_known_raw_neighbor_scores_across_alphabets(self):
        # Values transcribed from the inherited numeric tables. Direction is
        # preserved: M1(A,L)=1.575 whereas M1(L,A)=1.392.
        expected = {
            "AAAAA": 2.734,
            "PPPPP": 0.800,
            "CCCCC": 0.724,
            "GGGGG": 0.370 + 0.400 + 0.600,
            "LLLLL": 1.719 - 0.350 - 0.450,
            "YYYYY": 1.074 - 0.200 - 0.500,
            "ALALA": (1.575 + 1.392 + 1.575 + 1.392) / 4 + (-0.250 + 0.100) / 2,
            "AAAAL": (2.734 + 2.734 + 2.734 + 1.575) / 4 + (0.000 - 0.250) / 2 - 0.200,
        }
        for sequence, value in expected.items():
            minimal = list(dict.fromkeys(sequence))
            for alphabet in (FULL_ALPHABET, FULL_ALPHABET[::-1], minimal, minimal[::-1]):
                with self.subTest(sequence=sequence, alphabet=alphabet):
                    builder = make_builder(alphabet, len(sequence))
                    codes = [alphabet.index(aa) for aa in sequence]
                    self.assertAlmostEqual(builder._raw_helix_neigh(codes), value, places=12)

    def test_invalid_requested_alphabets_are_rejected(self):
        for alphabet in ([], ["A", "A"], ["A", "X"]):
            with self.subTest(alphabet=alphabet), self.assertRaises(ValueError):
                make_builder(alphabet)

    def test_matrix_change_invalidates_calibration_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            original = make_builder(FULL_ALPHABET)
            original.kwargs["stats_cache_dir"] = directory
            matrix, symbols = loader._load_third_neighbors_matrix_file()
            changed = matrix.copy()
            changed[0, 0] += 0.125
            with patch("core.hamiltonian_builder._load_third_neighbors_matrix_file", return_value=(changed, symbols)):
                updated = make_builder(FULL_ALPHABET)
            updated.kwargs["stats_cache_dir"] = directory
            self.assertNotEqual(original.matrix_hashes["M3"], updated.matrix_hashes["M3"])
            self.assertNotEqual(original._stats_cache_path(100), updated._stats_cache_path(100))
            self.assertGreater(original.stats_cache_version, 3)

    def test_fauchere_pliska_scale_and_descriptor_agree(self):
        from app.engine import HYDROPHOBICITY_FP

        builder = make_builder(FULL_ALPHABET)
        self.assertEqual(builder.hydro["S"], -0.040)
        self.assertEqual(builder.hydro["T"], 0.260)
        self.assertEqual(builder.hydro, HYDROPHOBICITY_FP)


class MatrixFormatTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.path = Path(self.temporary.name) / "matrix.txt"
        self.matrix = np.arange(400, dtype=float).reshape(20, 20)

    def write_matrix(self, symbols=None, rows=None):
        symbols = FULL_ALPHABET if symbols is None else symbols
        rows = self.matrix.tolist() if rows is None else rows
        self.path.write_text(
            " ".join(symbols) + "\n" + "\n".join(" ".join(map(str, row)) for row in rows) + "\n",
            encoding="utf-8",
        )

    def load_matrix(self):
        with patch.object(loader, "_construct_resource_path", return_value=str(self.path)):
            return loader._load_labeled_matrix_file("matrix.txt")

    def test_valid_matrix_round_trip(self):
        self.write_matrix()
        matrix, symbols = self.load_matrix()
        np.testing.assert_array_equal(matrix, self.matrix)
        self.assertEqual(symbols, FULL_ALPHABET)

    def test_invalid_headers_are_rejected(self):
        headers = (
            FULL_ALPHABET[:-1],
            FULL_ALPHABET + ["Y", "V"],
            ["A"] + FULL_ALPHABET[:-1],
            ["X"] + FULL_ALPHABET[1:],
        )
        for symbols in headers:
            with self.subTest(symbols=symbols):
                self.write_matrix(symbols=symbols)
                with self.assertRaises(ValueError):
                    self.load_matrix()

    def test_missing_extra_and_ragged_rows_are_rejected(self):
        rows = self.matrix.tolist()
        for malformed in (rows[:-1], rows + [rows[0]], [rows[0][:-1]] + rows[1:], [rows[0] + [0]] + rows[1:]):
            self.write_matrix(rows=malformed)
            with self.assertRaises(ValueError):
                self.load_matrix()

    def test_nonfinite_and_nonnumeric_values_are_rejected(self):
        for value in ("nan", "inf", "-inf", "bad"):
            with self.subTest(value=value):
                rows = self.matrix.tolist()
                rows[4][7] = value
                self.write_matrix(rows=rows)
                with self.assertRaises(ValueError):
                    self.load_matrix()


if __name__ == "__main__":
    unittest.main()
