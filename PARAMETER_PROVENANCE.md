# Parameter provenance and candidate corrections

**Publication candidate: the M3/M4 header reconstruction below is provisional.**
It is an inference from the inherited files, not an independently validated
assignment from a primary scientific source. Author confirmation of the original
labelled table/image and its bibliographic source is required before accepting
this candidate as a validated parameter release.

## Inherited matrices

The numerical tables were inherited from upstream commit
`61e9aadf8bfc960f9474a22f6aa63090f2fbded3`. The research workspace and its JCTC
validation snapshot contain the same malformed M3/M4 headers; agreement between
those copies is not independent evidence of scientific provenance.

| File | Header and interpretation in this candidate |
| --- | --- |
| `mj_matrix.txt` | Original 20 unique residue labels retained; original upper-triangle symmetrization retained. |
| `helix_pairs_prop.txt` (M1) | Original 20 unique residue labels retained; directional entries retained. The inherited use at separation 1 and positive sign are unchanged. |
| `helix_sd_i_3.txt` (M3) | Provisional 20-label reconstruction; directional entries retained. |
| `helix_sd_i_4.txt` (M4) | Provisional 20-label reconstruction; directional entries retained. |

The inherited M3/M4 header contained 22 labels for a 20 x 20 matrix:

```text
A G S T N D Q E H K R Y F W V I M L Y V C P
```

Keeping the first occurrence of each label gives this candidate order:

```text
A G S T N D Q E H K R Y F W V I M L C P
```

This removes only the repeated Y and V after L. It gives all 20 standard amino
acids, preserves the first 18 unique labels, and assigns C/P to the two final
rows/columns, which are all zero. Exchanging those two final labels would not
change numerical results. These facts support the reconstruction but cannot
prove the intended mapping or the numerical parameter source. The inherited
generation scripts mention an original image without identifying that image.

No numerical table values were edited. The complete bytes after the first
newline in each file match the upstream snapshot. Their SHA-256 digests are:

```text
mj_matrix.txt         8b988ae7a7faa0ead7c3c7542da5a5ea3006bc5a2a6b1593037ee7261fac31f6
helix_pairs_prop.txt  42bb5b35be6e6947f22c6e1da1b29b8fbea9597e53783bff639c1cb219de4b31
helix_sd_i_3.txt      d2dc9f7b8a296036fda5968f80efce68d28b5c354db02ca6b6a27f6b25ae6aa5
helix_sd_i_4.txt      a6bca28056910c9e8f0431878e19b0d9970b34e9c86c873725fbdcda9438c46f
```

## Confirmed hydrophobicity correction

The inherited implementation transposed the Fauchère–Pliska values for serine
and threonine. This candidate uses S = -0.040 and T = 0.260 in the scoring engine
and sequence descriptors, matching the scale published in the official
[HELIQUEST methods documentation](https://heliquest.ipmc.cnrs.fr/HelpProcedure.htm)
(checked 29 September 2026). This correction is independently supported by that
source and is separate from the provisional M3/M4 reconstruction.

## Indexing, cache invalidation, and limits of the checks

All four matrices are now selected by residue label on both axes for every
requested alphabet. Previously only MJ was remapped, so M1/M3/M4 scores could
change with the order or subset of the same chemical alphabet. The loader now
rejects duplicate/unknown labels, dimensions other than 20 x 20, and nonfinite
values. It does not silently repair future input files.

Calibration cache version 4 includes the labelled numerical matrix digests in
addition to the existing physicochemical tables. Digests cover the full loaded
matrix, with MJ hashed after its inherited symmetrization. Earlier calibration
caches cannot be reused for this candidate. Scores and rankings can change from
the earlier implementation; old and candidate outputs should not be mixed.

The regression tests cover full, permuted, and subset alphabets, directional
entries, explicit raw-score examples, malformed input rejection, scale
consistency, and cache invalidation. They verify implementation consistency,
not the provenance or predictive validity of the inherited matrices. Finite
Monte Carlo calibration still samples integer residue codes; permuting the
alphabet can change the sampled chemical decoys at a fixed seed, even when raw
scores are identical. No sign, interaction separation, or normalization formula
was changed; scientific interpretation of the inherited M1 term remains an
author review item alongside confirmation of the M3/M4 source.
