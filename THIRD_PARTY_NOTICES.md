# Third-party notices and scientific sources

The root MIT license applies to original project code and documentation. It does not replace notices or licenses on included third-party material.

## Included software

- **3Dmol.js 2.4.2**, copyright University of Pittsburgh and contributors, is distributed under BSD-3-Clause. Its license includes notices for GLmol, Three.js and jQuery. The complete upstream text is preserved in [LICENSES/3Dmol.txt](LICENSES/3Dmol.txt), with the companion bundle notice in `app/static/3Dmol-min.js.LICENSE.txt`. The bundle also contains pako (MIT and zlib), UPNG.js (MIT), iobuffer (MIT) and netcdfjs (MIT); their notices are reproduced in [LICENSES/3Dmol-bundled-dependencies.txt](LICENSES/3Dmol-bundled-dependencies.txt).
- `Code/data_loaders/energy_matrix_loader.py` derives from the IBM Qiskit protein-folding data loaders and carries the original **IBM 2021, 2022 / Apache-2.0** notice together with a modification notice. `Code/core/resources/mj_matrix.txt` is the Miyazawa–Jernigan table distributed with Qiskit under the same license. `Code/data_loaders/__init__.py` and `Code/core/resources/__init__.py` retain their original IBM notices. The Apache license is included in [LICENSES/Apache-2.0.txt](LICENSES/Apache-2.0.txt).
- Python dependencies are installed separately from the package index and retain their respective upstream licenses; they are not vendored in this repository.

## Parameter provenance

References for the scoring components are listed below. All tables are looked up by residue label; their headers are validated when they are loaded. This validation checks file structure, not agreement with a primary numerical source. The source and preprocessing of M1/M3/M4 remain unverified, and the M3/M4 residue ordering is provisional; see the [model limitations](README.md#model-scope-and-limitations).

- `mj_matrix.txt`: Miyazawa S, Jernigan RL. *Residue–residue potentials with a favorable contact pair term and an unfavorable high packing density term, for simulation and threading*. J Mol Biol **256**, 623–644 (1996), Table 3. [doi:10.1006/jmbi.1996.0114](https://doi.org/10.1006/jmbi.1996.0114).
- `helix_pairs_prop.txt`: supplied as i,i+1 residue-pair propensities. The implementation adds its entries directly to the M3/M4 contributions, then standardizes the combined neighbor term. The relationship between these numerical values and the normalized M⁽¹⁾ matrix of Eq. 3 has not been independently verified. A partial check against the FH counts and AH frequencies in Nacar's supplement reconstructs all 20 entries in the first row, but several checked entries differ: source-derived/repository values are R,M 1.244/1.204, N,A 0.817/0.882, D,A 0.721/1.271 and P,P 0/0.800. Source reconstruction uses `round(FH[j,i] * 1350293² / (306352 * AH[j] * AH[i]), 3)`; the source labels rows j and columns i. This does not establish the full table's derivation, sequence direction or its use at separation i,i+1. The stored values are retained until those choices are documented. Primary reference: Nacar C. *Propensities of Amino Acid Pairings in Secondary Structure of Globular Proteins*. Protein J **39**, 21–32 (2020). [doi:10.1007/s10930-020-09880-6](https://doi.org/10.1007/s10930-020-09880-6), [primary supplement](https://static-content.springer.com/esm/art%3A10.1007%2Fs10930-020-09880-6/MediaObjects/10930_2020_9880_MOESM1_ESM.pdf).
- `helix_sd_i_3.txt` and `helix_sd_i_4.txt`: supplied as i,i+3 and i,i+4 side-chain interaction energies (rows: residue i, columns: residue i+k; negative entries reduce the raw neighbor term). The provisional order assigns the zero final two rows and columns to Cys and Pro.
- Helix propensities in the scoring code: Pace CN, Scholtz JM. *A helix propensity scale based on experimental studies of peptides and proteins*. Biophys J **75**, 422–427 (1998).
- Hydrophobicity values in the scoring code: Fauchère JL, Pliska V. *Hydrophobic parameters II of amino acid side-chains from the partitioning of N-acetyl-amino acid amides*. Eur J Med Chem **18** (1983).

Named peptide examples in the interface are demonstration inputs.

## Branding

The SIMBIOS and Universidade de Santiago de Compostela names and logos identify the project and institution. They are excluded from the MIT grant; this distribution grants no trademark rights or endorsement. Replace these assets and institutional presentation when deploying an unrelated service.
