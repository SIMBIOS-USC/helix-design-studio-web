#!/usr/bin/env python3
"""Reconstruct Nacar's helical pairing propensities from published counts.

Supplement 1 supplies hydrogen-bonded pair counts F_H[j,i] and amino-acid
frequencies A_H. Methods section 2.4 defines P_o = F_H[j,i]/N_pairs and
P_r = A_H[j]*A_H[i]/(N*(N-1)); Figure 3 reports P_o/P_r to three decimals.
Source rows identify the later donor j and columns the earlier acceptor i.
The runtime table transposes these axes to use earlier residues as rows.

Only the Python standard library is required. This operation verifies the
tabulated coefficients; their application at separation 1 and negative sign
in the model are separate conventions.
"""

import argparse
import csv
from decimal import Decimal, ROUND_HALF_UP, localcontext
import hashlib
import json
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
SOURCE_DIRECTORY = REPOSITORY / "tools" / "data" / "parameter_sources"
RESOURCES_DIRECTORY = REPOSITORY / "Code" / "core" / "resources"
RESIDUES = set("ARNDCQEGHILKMFPSTWYV")
N_RESIDUES = 1_350_293
N_PAIRS = 306_352
PRINTED_PRECISION = Decimal("0.001")
SOURCE = {
    "doi": "10.1007/s10930-020-09880-6",
    "counts_location": "Supplement 1: helical F_H[j,i] and A_H counts",
    "formula_location": "Methods section 2.4",
    "propensities_location": "Figure 3, printed page 25",
    "supplement_url": (
        "https://static-content.springer.com/esm/art%3A10.1007%2Fs10930-020-09880-6/"
        "MediaObjects/10930_2020_9880_MOESM1_ESM.pdf"
    ),
    "supplement_sha256": "8ee7c96cf91fa459677d912abcccfba6361cba818f1c72fdb7ded34d1018c50c",
    "source_axes": {"rows": "later donor j", "columns": "earlier acceptor i"},
    "runtime_axes": {"rows": "earlier residue i", "columns": "later residue i+k"},
    "source_separation": 4,
    "runtime_separation": 1,
    "mapping": "transpose",
    "formula": "(F_H[j,i]/N_pairs)/(A_H[j]*A_H[i]/(N*(N-1)))",
    "rounding": "three decimal places, Decimal ROUND_HALF_UP",
}


def _labels(labels):
    if len(labels) != 20 or len(set(labels)) != 20 or set(labels) != RESIDUES:
        raise ValueError("Each axis must contain the 20 unique standard residue codes")
    return labels


def _integer(token):
    value = int(token)
    if value < 0:
        raise ValueError("Counts must be non-negative integers")
    return value


def _csv_matrix(path, convert):
    with Path(path).open(newline="", encoding="utf-8") as handle:
        records = list(csv.reader(handle))
    if len(records) != 21:
        raise ValueError("The source table must have one header and 20 rows")
    columns = _labels(records[0][1:])
    rows = _labels([record[0] for record in records[1:]])
    values = {}
    for row, record in zip(rows, records[1:]):
        if len(record) != 21:
            raise ValueError("Each source row must contain 20 values")
        for column, token in zip(columns, record[1:]):
            values[row, column] = convert(token)
    return rows, columns, values


def load_counts(source_directory=SOURCE_DIRECTORY):
    directory = Path(source_directory)
    rows, columns, pairs = _csv_matrix(directory / "nacar_fh_counts.csv", _integer)
    with (directory / "nacar_residue_counts.csv").open(newline="", encoding="utf-8") as handle:
        records = list(csv.reader(handle))
    if records[0] != ["residue", "helix_count"] or any(len(record) != 2 for record in records[1:]):
        raise ValueError("Residue count CSV requires residue and helix_count columns")
    _labels([record[0] for record in records[1:]])
    frequencies = {residue: _integer(count) for residue, count in records[1:]}
    if any(count == 0 for count in frequencies.values()):
        raise ValueError("All amino-acid frequencies must be positive")
    if sum(pairs.values()) != N_PAIRS or sum(frequencies.values()) != N_RESIDUES:
        raise ValueError("Transcribed count totals differ from the published totals")
    return rows, columns, pairs, frequencies


def reconstruct(pairs, frequencies):
    """Return the printed-precision propensities on the published j,i axes."""
    with localcontext() as context:
        context.prec = 40
        return {
            (row, column): (
                Decimal(count) * Decimal(N_RESIDUES * (N_RESIDUES - 1))
                / Decimal(N_PAIRS * frequencies[row] * frequencies[column])
            ).quantize(PRINTED_PRECISION, rounding=ROUND_HALF_UP)
            for (row, column), count in pairs.items()
        }


def _runtime_matrix(path):
    lines = Path(path).read_text(encoding="utf-8").splitlines()
    header = lines[0].split()
    labels = _labels(list(header[0]) if len(header) == 1 else header)
    if len(lines) != 21:
        raise ValueError("The runtime table must have one header and 20 rows")
    values = {}
    for row, line in zip(labels, lines[1:]):
        tokens = line.split()
        if len(tokens) != 20:
            raise ValueError("Each runtime row must contain 20 values")
        for column, token in zip(labels, tokens):
            values[row, column] = Decimal(token)
    return values


def _differences(expected, actual, row_label, column_label):
    return [{row_label: row, column_label: column,
             "reconstructed": str(value), "table_value": str(actual[row, column])}
            for (row, column), value in expected.items() if actual[row, column] != value]


def verify(source_directory=SOURCE_DIRECTORY, resources_directory=RESOURCES_DIRECTORY):
    directory = Path(source_directory)
    rows, columns, pairs, frequencies = load_counts(directory)
    reconstructed = reconstruct(pairs, frequencies)
    _, _, printed = _csv_matrix(directory / "nacar2020_figure3.csv", Decimal)
    runtime_path = Path(resources_directory) / "helix_pairs_prop.txt"
    runtime = _runtime_matrix(runtime_path)
    printed_differences = _differences(reconstructed, printed, "source_row_j", "source_column_i")
    transposed = {(column, row): value for (row, column), value in reconstructed.items()}
    runtime_differences = _differences(transposed, runtime, "residue_i", "residue_i_plus_k")
    paths = [directory / "nacar_fh_counts.csv", directory / "nacar_residue_counts.csv",
             directory / "nacar2020_figure3.csv", runtime_path]
    return {
        "source": SOURCE,
        "source_row_order": "".join(rows), "source_column_order": "".join(columns),
        "pair_count_total": sum(pairs.values()), "residue_count_total": sum(frequencies.values()),
        "reconstructed_cells": len(reconstructed),
        "published_cells_matching": 400 - len(printed_differences),
        "published_cell_differences": printed_differences,
        "runtime_cells_matching": 400 - len(runtime_differences),
        "runtime_cell_differences": runtime_differences,
        "input_sha256": {path.name: hashlib.sha256(path.read_bytes()).hexdigest() for path in paths},
        "verified": not printed_differences and not runtime_differences,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, default=SOURCE_DIRECTORY)
    parser.add_argument("--resources-dir", type=Path, default=RESOURCES_DIRECTORY)
    parser.add_argument("--output", type=Path, help="Save the complete JSON comparison")
    args = parser.parse_args()
    result = verify(args.source_dir, args.resources_dir)
    serialized = json.dumps(result, indent=2, ensure_ascii=False) + "\n"
    if args.output:
        args.output.write_text(serialized, encoding="utf-8")
    print(serialized, end="")
    return 0 if result["verified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
