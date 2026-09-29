# Third-party notices and scientific sources

The root MIT license applies to original project code and documentation. It does not replace notices or licenses on included third-party material.

## Included software

- **3Dmol.js 2.4.2**, copyright University of Pittsburgh and contributors, is distributed under BSD-3-Clause. Its license includes notices for GLmol, Three.js and jQuery. The complete upstream text is preserved in [LICENSES/3Dmol.txt](LICENSES/3Dmol.txt), with the companion bundle notice in `app/static/3Dmol-min.js.LICENSE.txt`.
- `Code/data_loaders/__init__.py` and `Code/core/resources/__init__.py` retain their original **IBM 2021, 2022 / Apache-2.0** notices. The Apache license is included in [LICENSES/Apache-2.0.txt](LICENSES/Apache-2.0.txt).
- Python dependencies are installed separately from the package index and retain their respective upstream licenses; they are not vendored in this repository.

## Parameter provenance

References for the scoring components are listed below. The numerical provenance of the helix-neighbor tables is not independently verified; see [KNOWN_ISSUES.md](KNOWN_ISSUES.md).

- `mj_matrix.txt`: Miyazawa S, Jernigan RL. *Residue–residue potentials with a favorable contact pair term and an unfavorable high packing density term, for simulation and threading*. J Mol Biol **256**, 623–644 (1996), Table 3. [doi:10.1006/jmbi.1996.0114](https://doi.org/10.1006/jmbi.1996.0114).
- Helix-neighbor tables: related reference: Nacar C. *Propensities of Amino Acid Pairings in Secondary Structure of Globular Proteins*. Protein J **39**, 21–32 (2020). [doi:10.1007/s10930-020-09880-6](https://doi.org/10.1007/s10930-020-09880-6). The exact transformation and row/column provenance of all three supplied tables remain to be verified; see `KNOWN_ISSUES.md`.
- Helix propensities in the scoring code: Pace CN, Scholtz JM. *A helix propensity scale based on experimental studies of peptides and proteins*. Biophys J **75**, 422–427 (1998).
- Hydrophobicity values in the scoring code: Fauchère JL, Pliska V. *Hydrophobic parameters II of amino acid side-chains from the partitioning of N-acetyl-amino acid amides*. Eur J Med Chem **18** (1983).

Named peptide examples in the interface are demonstration inputs.

## Branding

The SIMBIOS and Universidade de Santiago de Compostela names and logos identify the project and institution. They are excluded from the MIT grant; this distribution grants no trademark rights or endorsement. Replace these assets and institutional presentation when deploying an unrelated service.
