# Source provenance

This web-only source distribution was assembled on 2026-09-29 from the project maintainers' source snapshot `61e9aadf8bfc960f9474a22f6aa63090f2fbded3` (source commit dated 2026-08-26). It starts a new repository history and is not a copy of the research workspace or of private repository history.

## Initial distribution changes (4cbc8f5)

- Selected the API, engine, static interface, required assets and four parameter tables explicitly.
- Removed unused Pauli-operator/quantum Hamiltonian construction methods and their helper; retained classical scoring/calibration routines. Removed duplicate loader definitions.
- Kept all four numeric resource files unchanged. See `KNOWN_ISSUES.md` for the unresolved matrix-indexing defect.
- Restricted automatic runtime discovery to the bundled `Code/` directory; an explicit `QFOLD_CODE_DIR` override remains available.
- Made calibration and optional usage-log paths configurable. Disabled application usage logging and persistent visitor identifiers by default; removed raw internal error text from public health diagnostics.
- Updated the interface's method citation from bioRxiv to the published JCTC article.
- Added installation instructions, pinned runtime dependencies, a portable container definition, tests, citation metadata, MIT licensing for original code and third-party notices.

No production deployment was performed. The public website did not provide its deployed source revision during verification, so equivalence with its running backend is not asserted.

## Molecular viewer

`app/static/3Dmol-min.js` matches the official **3Dmol.js 2.4.2** npm distribution byte for byte:

```text
SHA256 06b6d2fc7d418e8bef62a32cca44373ff21b2b787bd377613dd4039394a4e9ff
size   513564 bytes
```

The bundle, its companion notice and full license were checked against the pinned package at `https://cdn.jsdelivr.net/npm/3dmol@2.4.2/`. The companion notice and full license are included in this repository. Runtime requests to that CDN are unnecessary.

## Review candidate 0.2.0rc1

The candidate maps all matrices by residue identity and rejects malformed tables. The two duplicated Y/V header entries in M3/M4 were removed provisionally, without changing any numeric entries; source confirmation remains pending and is documented in PARAMETER_PROVENANCE.md. Calibration caches now identify parameter content. The alphabet-order regression is expected to pass.

Specificity and cross-environment result preparation now reuse the actual search builders, avoiding post-search changes of calibration and orientation. Responses record calibration and reference seeds and orientation mode. These corrections change scores relative to the initial snapshot and require a fresh evaluation. The production website has not been updated by this repository preparation.

The Ser/Thr entries in the Fauchère–Pliska hydrophobicity scale were corrected to S=-0.040 and T=0.260 in both scoring and descriptor code, matching the official HELIQUEST parameter table. This is separate from the provisional M3/M4 header reconstruction.
