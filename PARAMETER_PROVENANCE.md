# Parameter definitions and provenance

The four parameter files are in `Code/core/resources/`. Rows and columns use the residue labels in each file's header.

| File | Computational interpretation |
| --- | --- |
| `mj_matrix.txt` | Miyazawa–Jernigan pair matrix; the upper triangle is symmetrized and the diagonal retained. |
| `helix_pairs_prop.txt` (M1) | Directional entries at sequence separation 1, added with a positive sign. |
| `helix_sd_i_3.txt` (M3) | Directional entries at separation 3; provisional residue order. |
| `helix_sd_i_4.txt` (M4) | Directional entries at separation 4; provisional residue order. |

Scientific references are listed in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). The exact numerical source and transformation of the helix-neighbor tables, including the sign convention for M1, have not been independently verified.

## M3/M4 residue order

Both 20 × 20 matrices use this provisional row/column order:

```text
A G S T N D Q E H K R Y F W V I M L C P
```

The assignment is inferred from a 22-label header containing duplicate Y and V entries, retaining the first occurrence of each label. The final two rows and columns are zero, so exchanging C and P has no numerical effect. This internal consistency does not establish the intended scientific assignment. Results are conditional on the stated ordering; calculations must be repeated if it changes.

## Table checksums

SHA-256 identifiers for the numerical rows, excluding the first header line:

```text
mj_matrix.txt         8b988ae7a7faa0ead7c3c7542da5a5ea3006bc5a2a6b1593037ee7261fac31f6
helix_pairs_prop.txt  42bb5b35be6e6947f22c6e1da1b29b8fbea9597e53783bff639c1cb219de4b31
helix_sd_i_3.txt      d2dc9f7b8a296036fda5968f80efce68d28b5c354db02ca6b6a27f6b25ae6aa5
helix_sd_i_4.txt      a6bca28056910c9e8f0431878e19b0d9970b34e9c86c873725fbdcda9438c46f
```

## Hydrophobicity

Scoring and sequence descriptors use the Fauchère–Pliska scale, with Ser = −0.040 and Thr = 0.260. The scale values are listed in the [HELIQUEST methods documentation](https://heliquest.ipmc.cnrs.fr/HelpProcedure.htm).

## Loading, calibration and tests

All four matrices are indexed by residue identity for the requested alphabet. The loader requires 20 unique standard residue labels and a finite 20 × 20 numeric array. MJ is symmetrized; the helix-neighbor matrices retain their directionality.

Calibration-cache version 4 includes hashes of the loaded matrices, physicochemical tables and request settings. MJ is hashed after symmetrization. Reproduce calculations using matching source, parameters and calibration settings.

Tests cover full, permuted and subset alphabets, directional entries, raw-score examples, invalid inputs and cache invalidation. These tests establish implementation consistency rather than scientific validity. At a fixed seed, changing alphabet order can change the sampled chemical decoys, even when raw scores are invariant.
