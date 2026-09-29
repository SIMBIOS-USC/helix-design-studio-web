"""Load residue-labelled parameter matrices with strict format checks."""

import os
from typing import List, Tuple

import numpy as np


_STANDARD_RESIDUES = frozenset("ACDEFGHIKLMNPQRSTVWY")


def _construct_resource_path(filename: str) -> str:
    return os.path.realpath(
        os.path.join(os.path.dirname(__file__), "..", "core", "resources", filename)
    )


def _load_labeled_matrix_file(filename: str) -> Tuple[np.ndarray, List[str]]:
    """Read exactly one finite 20 x 20 matrix, with 20 unique residue labels.

    Rows and columns use the same header order. Blank lines after the header
    are allowed; missing, extra, or ragged data rows are rejected.
    """
    path = _construct_resource_path(filename)
    with open(path, "r", encoding="utf-8") as handle:
        symbols = handle.readline().split()
        if len(symbols) != 20:
            raise ValueError(f"{filename}: expected 20 residue labels, got {len(symbols)}")
        if len(set(symbols)) != len(symbols):
            raise ValueError(f"{filename}: residue labels must be unique")
        if set(symbols) != _STANDARD_RESIDUES:
            raise ValueError(f"{filename}: labels must contain the 20 standard amino acids")

        rows = [line.split() for line in handle if line.strip()]

    if len(rows) != 20 or any(len(row) != 20 for row in rows):
        raise ValueError(f"{filename}: expected exactly 20 data rows with 20 values each")
    try:
        matrix = np.array(rows, dtype=float)
    except ValueError as error:
        raise ValueError(f"{filename}: matrix values must be numeric") from error
    if not np.isfinite(matrix).all():
        raise ValueError(f"{filename}: matrix values must be finite")
    return matrix, symbols


def _load_energy_matrix_file() -> Tuple[np.ndarray, List[str]]:
    """Load MJ and symmetrize its upper triangle, preserving the diagonal."""
    raw, symbols = _load_labeled_matrix_file("mj_matrix.txt")
    upper = np.triu(raw, k=1)
    return np.diag(np.diag(raw)) + upper + upper.T, symbols


def _load_first_neighbors_matrix_file() -> Tuple[np.ndarray, List[str]]:
    """Load the directional helix-pair propensity matrix."""
    return _load_labeled_matrix_file("helix_pairs_prop.txt")


def _load_third_neighbors_matrix_file() -> Tuple[np.ndarray, List[str]]:
    """Load i,i+3 parameters; see PARAMETER_PROVENANCE.md for header status."""
    return _load_labeled_matrix_file("helix_sd_i_3.txt")


def _load_fourth_neighbors_matrix_file() -> Tuple[np.ndarray, List[str]]:
    """Load i,i+4 parameters; see PARAMETER_PROVENANCE.md for header status."""
    return _load_labeled_matrix_file("helix_sd_i_4.txt")
