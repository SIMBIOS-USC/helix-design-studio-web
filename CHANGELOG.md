# Changelog

Versions are tagged in the repository (`v0.1.0`, `v1.0.0`). Each result reports the version, the parameter-set checksum, the complete request and the seed protocol under `provenance`; `GET /api/health` reports the running version. Results from different versions are not interchangeable.

## 1.0.0 — 2026-09-30

### Scoring model

Measured with `tools/measure_changes.py`: each row restores one 0.1.0 behaviour on top of 1.0.0, so its effect is isolated. Benchmark: the 12 reference peptides of the interface in the 5 environment presets (60 Z-scores, full alphabet, 1000 calibration decoys, 500 random references, seed 42); the mean percentile of 150 random 18-residue sequences in the anionic interface (50 for an unbiased score); and 6 designs of 18 residues with fixed seeds. For comparison, designs from the same model with other seeds differ in 27% of positions.

| Change | 0.1.0 | 1.0.0 | Reason | Median / max \|ΔZ\| | Spearman of Z | Design positions changed |
|---|---|---|---|---|---|---|
| Helix-neighbor tables | M1 (i,i+1), M3 (i,i+3) and M4 (i,i+4) indexed by position in the requested alphabet | Indexed by residue label, as the Miyazawa–Jernigan table already was | Positional indexing assigned the parameters of other residues whenever the alphabet order differed from the file order (always for M3/M4) | 0.111 / 0.284 | 0.993 | 56% |
| i,i+1 pair term | Propensity *P* added to the minimized energy | Energy −ln *P* | *P* > 1 marks pairs over-represented in helices; adding *P* favoured helix-breaking pairs (in this term, poly-Gly scored better than poly-Ala) | 0.375 / 3.349 | 0.802 | 81% |
| Fauchère–Pliska hydrophobicity | Ser 0.26, Thr −0.04 | Ser −0.04, Thr 0.26 | Values were swapped relative to the published scale | 0.033 / 0.142 | 0.999 | 31% |
| Pace–Scholtz helix propensity | Pro 3.15 | Pro 3.16 | Published value | 0.001 / 0.002 | 1.000 | 3% |
| Miyazawa–Jernigan contacts | i to i+1 (d_max = 1) | i to i+1…i+4 (d_max = 4) | Contact range of the method article | 0.001 / 0.007 | 1.000 | 16% |
| Orientation in scoring | Only the query was oriented by its hydrophobic moment; calibration decoys and random references were evaluated at the query's phase | Query, calibration decoys and random references are oriented by the same rule | The query had an orientation advantage over its background: random sequences scored at a mean percentile of 25.9 instead of 50 (51.7 in 1.0.0) | 0.328 / 1.592 | 0.925 | scoring only |

All changes together (complete 0.1.0 source tree against 1.0.0, same benchmark): median |ΔZ| 0.632, maximum 2.817 (Poly-Ala 15, homogeneous apolar), Spearman 0.818; mean percentile of random sequences 22.1 in 0.1.0 and 51.7 in 1.0.0. These figures compare software versions on this benchmark; assessing results obtained with earlier versions requires their original inputs and settings.

Also in the model:

- The i,i+3 and i,i+4 table headers listed 22 labels for 20 columns (Y and V repeated after L). They now list 20 labels: the first 18 as given, and C and P for the last row and column, which are all zero.
- Searches (Design, Family, Cross-Design, Specificity) optimize the fixed-geometry Hamiltonian, as in 0.1.0. Design results now report both the design-geometry energy and the standard score. Cross-Design and Specificity previously reported energies and margins from a separate re-scoring (other calibration seeds, sample sizes and orientation); they now report the optimized quantities, and rankings use them.
- The penetration scan uses the same seed at every depth (previously one seed per depth).
- Calibration, random references and searches draw from independent random streams (`SeedSequence`); 0.1.0 used seed offsets that could coincide. Results for a given seed therefore differ from 0.1.0 within sampling noise.
- The calibration no longer evaluates a hydrophobic-moment term that was not part of the energy (no effect on results).
- Calibration caches are keyed by the parameter tables; caches from 0.1.0 are not reused.

### API and interface

- Stricter validation: standard residues only, the sequence must use the allowed residues, known presets and membrane charges, bounded angles, weights and seeds, no NaN/Infinity, no unknown fields, at most four off-target environments. Physically identical target and off-target environments are rejected. `wheel_phase_deg` is no longer accepted. `/api/design` uses the 20-residue alphabet by default, like the interface.
- Client errors return 400/422 and internal failures 500; `/api/health` only reports server-side failures and includes the version.
- Calculations have a concurrency limit, a time limit and a request-size limit, and stop when the client disconnects (including streams and calibration).
- Results include `provenance` and `calibration` blocks; the interface shows the settings used and exports them in JSON downloads.
- Interface fixes: cancel buttons, partial results kept on stop or error, consistent meaning of the penetration slider, helical wheel drawn with the model's angular sector, escaped error messages, accessibility and layout fixes. Two reference peptides whose sequences did not match the cited proteins (GB1 helix fragment, α3D segment) were removed and KALP17 follows the GKK(LA)ₙLKKA pattern.

### Reproducing the measurements

```bash
git worktree add ../helix-0.1.0 v0.1.0
python tools/measure_changes.py --baseline ../helix-0.1.0
```

## 0.1.0

Web application source before the changes above (API version 0.1.0, commit `a8cdf43`).
