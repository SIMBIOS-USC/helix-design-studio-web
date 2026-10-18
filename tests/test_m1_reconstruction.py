"""Independent reconstruction of published M1 coefficients from integer counts."""

import csv
from decimal import Decimal
from pathlib import Path
import shutil
import tempfile
import unittest

from tools import reconstruct_m1_from_counts as reconstruction


class M1ReconstructionTests(unittest.TestCase):
    def test_counts_reconstruct_source_and_runtime_at_printed_precision(self):
        result = reconstruction.verify()
        self.assertTrue(result["verified"], result)
        self.assertEqual(result["pair_count_total"], 306352)
        self.assertEqual(result["residue_count_total"], 1350293)
        self.assertEqual(result["published_cells_matching"], 400)
        self.assertEqual(result["runtime_cells_matching"], 400)

    def test_direction_and_zero_donor_proline_are_preserved(self):
        _, _, counts, frequencies = reconstruction.load_counts()
        source = reconstruction.reconstruct(counts, frequencies)
        self.assertEqual(source["A", "R"], Decimal("1.648"))
        self.assertEqual(source["R", "A"], Decimal("1.724"))
        self.assertTrue(all(source["P", residue] == 0 for residue in reconstruction.RESIDUES))
        self.assertNotEqual(source["A", "P"], 0)

    def test_an_untransposed_runtime_is_detected(self):
        with tempfile.TemporaryDirectory(prefix="helix-m1-direction-") as temporary:
            directory = Path(temporary)
            rows, columns, source = reconstruction._csv_matrix(
                reconstruction.SOURCE_DIRECTORY / "nacar2020_figure3.csv", Decimal)
            text = " ".join(columns) + "\n"
            text += "\n".join(" ".join(str(source[row, column]) for column in columns)
                              for row in rows) + "\n"
            (directory / "helix_pairs_prop.txt").write_text(text, encoding="utf-8")
            result = reconstruction.verify(resources_directory=directory)
            self.assertFalse(result["verified"])
            self.assertEqual(result["published_cells_matching"], 400)
            self.assertLess(result["runtime_cells_matching"], 400)

    def test_pair_count_error_is_detected_even_if_the_total_is_preserved(self):
        with tempfile.TemporaryDirectory(prefix="helix-m1-count-") as temporary:
            directory = Path(temporary)
            for path in reconstruction.SOURCE_DIRECTORY.glob("nacar*.csv"):
                shutil.copyfile(path, directory / path.name)
            path = directory / "nacar_fh_counts.csv"
            with path.open(newline="", encoding="utf-8") as handle:
                records = list(csv.reader(handle))
            records[1][1] = str(int(records[1][1]) + 100)
            records[1][2] = str(int(records[1][2]) - 100)
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle).writerows(records)
            result = reconstruction.verify(source_directory=directory)
            self.assertFalse(result["verified"])
            self.assertEqual(result["pair_count_total"], 306352)
            self.assertEqual(len(result["published_cell_differences"]), 2)
            self.assertEqual(len(result["runtime_cell_differences"]), 2)

    def test_missing_or_duplicated_residue_counts_are_rejected(self):
        with tempfile.TemporaryDirectory(prefix="helix-m1-frequency-") as temporary:
            directory = Path(temporary)
            for path in reconstruction.SOURCE_DIRECTORY.glob("nacar*.csv"):
                shutil.copyfile(path, directory / path.name)
            path = directory / "nacar_residue_counts.csv"
            records = path.read_text(encoding="utf-8").splitlines()
            records[-1] = records[-2]
            path.write_text("\n".join(records) + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unique standard residue"):
                reconstruction.verify(source_directory=directory)


if __name__ == "__main__":
    unittest.main()
