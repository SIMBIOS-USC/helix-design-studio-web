# Changelog

Use `v0.1.0` for the comparison baseline and `v1.0.0rc1` for this release candidate. Each result reports the version, the parameter-set checksum, the complete request and the seed protocol under `provenance`; `GET /api/health` reports the running version. Results from different versions are not interchangeable.

## 1.0.0rc1 — 2026-10-12

Candidate reviewed on 2026-10-12. The following measurements apply to the parameter files identified in the README and to the protocols below.

### Scoring model

Measured with `tools/measure_changes.py`: each row restores one 0.1.0 behaviour on top of the candidate, with the other settings held fixed. Benchmark: the 12 reference peptides of the interface in the 5 environment presets (60 Z-scores, full alphabet, 1000 calibration decoys, 500 random references, seed 42); the mean percentile of 150 random 18-residue sequences in the anionic interface (50 is the population expectation for a continuous, exchangeable query/reference distribution); and 6 designs of 18 residues with fixed seeds. The "design positions changed" column describes those six runs. It does not separate sensitivity to model parameters from search variability. The script also reports a separate search control: it changes search entropy from [seed,2] to [seed+100,2], while holding the calibration and random-reference streams fixed. This control changed 41% of design positions over the six runs; it is a descriptive sensitivity check, not an attribution of parameter effects.

| Change | 0.1.0 | 1.0.0rc1 | Reason | Median / max \|ΔZ\| | Spearman of Z | Design positions changed |
|---|---|---|---|---|---|---|
| Helix-neighbor tables | M1 (i,i+1), M3 (i,i+3) and M4 (i,i+4) indexed by position in the requested alphabet | Indexed by residue label, as the Miyazawa–Jernigan table already was | Positional indexing assigned the parameters of other residues whenever the alphabet order differed from the file order (always for M3/M4) | 0.173 / 0.499 | 0.990 | 47% |
| Fauchère–Pliska hydrophobicity | Ser 0.26, Thr −0.04 | Ser −0.04, Thr 0.26 | Values were swapped relative to the published scale | 0.033 / 0.145 | 0.999 | 59% |
| Pace–Scholtz helix propensity | Pro 3.15 | Pro 3.16 | Published value | 0.001 / 0.002 | 1.000 | 31% |
| Miyazawa–Jernigan contacts | i to i+1 (d_max = 1) | i to i+1…i+4 (d_max = 4) | Contact range of the method article (default d_max = 4) | 0.001 / 0.008 | 1.000 | 26% |
| Orientation in scoring | Only the query was oriented by its hydrophobic moment; calibration decoys and random references were evaluated at the query's phase | Query, calibration decoys and random references are oriented by the same rule | The query had an orientation advantage over its background: random sequences scored at a mean percentile of 24.5 (51.2 in the candidate) | 0.363 / 1.626 | 0.951 | scoring only |

The i,i+1, i,i+3 and i,i+4 tables are summed directly, with their stored sign and no per-matrix transformation; 0.1.0 used the same convention. The combined neighbor term is standardized after summation. Eq. 3 motivates an additive term, but does not verify the origin, preprocessing or residue ordering of the supplied numerical tables. Those source limitations remain open.

All changes together (complete 0.1.0 source tree against the candidate, same benchmark): median |ΔZ| 0.491, maximum 1.543 (Poly-Ala 15, neutral interface), Spearman 0.944; mean percentile of random sequences 22.1 in 0.1.0 and 51.2 in the candidate. The orientation comparison above produces the largest median score change among the individual variants in this benchmark. These figures compare software versions on the declared inputs; assessing results obtained with earlier versions requires their original inputs and settings.

Also in the model:

- The i,i+3 and i,i+4 table headers listed 22 labels for 20 columns (Y and V repeated after L). They now list 20 labels: the first 18 as given, and C and P for the last row and column, which are all zero.
- Searches (Design, Family, Cross-Design, Specificity) optimize the fixed-geometry Hamiltonian, as in 0.1.0. Design results now report both the design-geometry energy and the standard score. Cross-Design and Specificity previously reported energies and margins from a separate re-scoring (other calibration seeds, sample sizes and orientation); they now report the optimized quantities, and rankings use them.
- The penetration scan uses the same seed at every depth (previously one seed per depth).
- Calibration, random references and searches draw from independent random streams (`SeedSequence`); 0.1.0 used seed offsets that could coincide. The same user seed therefore selects different random draws across versions.
- The calibration no longer evaluates a hydrophobic-moment term that was not part of the energy (no effect on results).
- Calibration caches are keyed by the parameter tables; caches from 0.1.0 are not reused.

### API and interface

- Stricter validation: standard residues only, the sequence must use the allowed residues, known presets and membrane charges, bounded angles, weights and seeds, no NaN/Infinity, no unknown fields, at most four off-target environments. Physically identical target and off-target environments are rejected. `wheel_phase_deg` is no longer accepted. `/api/design` uses the 20-residue alphabet by default, like the interface.
- Client errors return 400/422 and internal failures 500; `/api/health` only reports server-side failures and includes the version.
- Calculations have a concurrency limit, a time limit and a request-size limit, and stop when the client disconnects (including streams and calibration).
- Results include `provenance` and `calibration` blocks. Stream metadata and every fully evaluated candidate carry provenance, so interrupted Family and Specificity downloads preserve the complete request and version. Both workflows export their partial JSON results after stop or error.
- Interface fixes: cancel buttons, partial results kept on stop or error, consistent meaning of the penetration slider, helical wheel drawn with the model's angular sector, escaped error messages, accessibility and layout fixes. Two reference peptides whose sequences did not match the cited proteins (GB1 helix fragment, α3D segment) were removed and KALP17 follows the GKK(LA)ₙLKKA pattern.

### Reproducing the measurements

```bash
git worktree add ../helix-0.1.0 v0.1.0
python tools/measure_changes.py --baseline ../helix-0.1.0
```

Add `--output model_change_comparison.json` to save the full input protocol, scores, designs and numerical summary, including the control that changes only the search stream.

### Comparison with the method article's Table 1

`tools/data/article_table1.json` contains seven exact sequence strings and fourteen Z-scores transcribed from the identified author manuscript's Table 1, with its DOI, source hash and the original table rows. GWALP23 is printed with length 23 but contains 21 residues in that sequence string; the rerun evaluates the printed sequence and records both lengths. Terminal modifications are omitted, as in the sequence-only model.

```bash
python tools/reproduce_article_table1.py --output article_table1_comparison.json
```

The declared rerun uses the full 20-residue alphabet in its recorded order, default weights, 2,000 calibration sequences, 10,000 random references per environment, seed 42 and the candidate's independent random streams. It compares direct M1 entries with −ln(M1), holding every other setting and the random sequences fixed while recalibrating the combined neighbor term. The output includes all fourteen scores for each convention, error summaries, input and source-file hashes, and runtime versions. Agreement with the rounded reference values does not establish the provenance of M1/M3/M4 or reproduce the historical calculation protocol.

The [archived numerical output](tools/data/article_table1_comparison.json) gives a mean absolute Z-score difference of **0.444852** for direct M1 entries and **1.313673** for −ln(M1), over the fourteen transcribed reference values. The maximum absolute differences are 0.670831 and 1.666431, respectively. These are comparisons with the rounded values under this rerun protocol, not estimates of experimental accuracy.

## 0.1.0

Web application source before the changes above (API version 0.1.0, commit `a8cdf43`).
