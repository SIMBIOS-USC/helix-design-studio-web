# Known issues and model limitations

## Helix-neighbor matrix indexing

`HamiltonianBuilder._load_matrices_from_files` maps the Miyazawa–Jernigan matrix by residue label but indexes M1/M3/M4 by positions in the requested residue alphabet. Changing alphabet order can therefore change the raw helix-neighbor score of the **same chemical sequence**. The design and scoring workflows use different default alphabets, which can affect comparisons, rankings and generated sequences.

The headers of `helix_sd_i_3.txt` and `helix_sd_i_4.txt` contain 22 symbols, including repeated Y and V, for 20 × 20 numeric matrices. Their row/column identities require verification against the parameter source. The expected-failure regression test documents the indexing defect. A fixed alphabet order allows reproducing a calculation but does not resolve the defect.

## Interpretation

- Scores measure compatibility with an assumed alpha-helical state. They do not establish folding, binding, membrane insertion, biological activity or experimental stability.
- PDB export constructs an ideal helix; it does not predict or minimize a structure.
- The Penetration scan changes an angular exposure mask. Its coordinate is not a physical insertion depth or a free-energy reaction coordinate.
- Calibration and reference distributions depend on sequence length, alphabet, environment, sample sizes and seeds. Scores from different protocols are not directly interchangeable.
- Execution tests do not establish parameter provenance or predictive validity.
