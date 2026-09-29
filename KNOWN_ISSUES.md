# Known scientific limitations

## Helix-neighbor matrix indexing

The current `HamiltonianBuilder._load_matrices_from_files` uses residue labels to map the Miyazawa–Jernigan matrix, but discards the labels supplied with the three helix-neighbor tables (`M1`, `M3`, `M4`). Those tables are subsequently indexed by positions in the request's residue alphabet. Consequently, changing the alphabet order can change the raw helix-neighbor score of the **same chemical sequence**. The default design and scoring workflows also use different alphabets, so this can affect comparisons between workflows, rankings and selected sequences.

Additionally, the headers of `helix_sd_i_3.txt` and `helix_sd_i_4.txt` contain 22 symbols, including repeated Y and V, while their numeric matrices are 20×20. Correct row/column identities must be established from the parameter source before changing the mapping. The runtime tables and scoring equations have therefore been preserved in this initial source distribution rather than silently reconstructed.

The expected-failure regression test documents the required alphabet-order invariance. Fixing this issue requires verifying table provenance/order, mapping by residue identity, checking all alphabets and rerunning the affected scientific comparisons. Keeping one alphabet order fixed is useful for reproducing this snapshot but **does not resolve the underlying defect**.

## Interpretation

- Scores concern compatibility of sequences with an assumed alpha-helical state. They do not establish folding, binding affinity, membrane insertion, biological activity or experimental stability.
- PDB export builds an ideal alpha helix using PeptideBuilder; it is not a structure prediction or an energy-minimized model.
- Calibration and random-reference distributions depend on length, alphabet, environment, sample size and seed. Scores from distinct protocols should not be assumed directly comparable.
- Automated execution tests and agreement with the source snapshot establish implementation continuity only. They do not validate the parameter tables or the scientific conclusions.
