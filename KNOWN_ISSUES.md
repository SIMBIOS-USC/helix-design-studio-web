# Model limitations and scoring conventions

## Parameter assignments

The M3/M4 residue order is provisional and has not been verified against a primary parameter source. Scores and designs depend on this assignment. The numerical source and transformation of the helix-neighbor tables, including the M1 sign and interaction separation, also require verification. See [PARAMETER_PROVENANCE.md](PARAMETER_PROVENANCE.md) for the implemented conventions and table identifiers.

## Search and rescoring

Specificity and cross-design outputs use the calibration and fixed geometry of their search objectives. Score and Compare align interfacial sequences to their hydrophobic moments. Rescoring a generated sequence through those workflows therefore changes the evaluation protocol. Different residue alphabets define different random-reference distributions.

## Interpretation

- Scores measure compatibility with an assumed alpha-helical state. They do not establish folding, binding, membrane insertion, biological activity or experimental stability.
- PDB export constructs an ideal helix; it does not predict or minimize a structure.
- The 21-point Penetration scan changes an angular exposure mask, with homogeneous endpoints. Its coordinate is not a physical insertion depth or a free-energy reaction coordinate.
- Calibration and reference sample sizes, alphabet order, orientation and seeds affect results. Finite-sampling variation can affect comparisons.
- Implementation tests do not establish parameter provenance or predictive validity.
