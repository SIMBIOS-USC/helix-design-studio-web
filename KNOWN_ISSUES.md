# Scientific status of review candidate 0.2.0rc1

## Parameter provenance requires author confirmation

The candidate repairs the helix-neighbor indexing defect in the initial public snapshot: matrices are now mapped using residue labels, with dimensions, unique labels and finite values checked. Raw-score alphabet-order invariance is covered by regression tests.

The supplied M3/M4 numeric arrays have 20 rows and columns, but the original headers contained 22 symbols, repeating Y and V. This candidate removes the duplicate header entries while preserving first occurrence and all numeric entries. This is a provisional reconstruction of internal file intent, not independent verification against an original scientific source. See PARAMETER_PROVENANCE.md for evidence and limitations. Author confirmation of the source/order remains necessary before submission. Recompute results if that interpretation changes.

## Search and rescoring conventions

Specificity and cross-design outputs now retain the calibration and fixed geometry used during search; their displayed energies reconstruct the recorded objective. Score and Compare automatically align interfacial sequences to their hydrophobic moments. Rescoring a generated sequence through those workflows therefore changes the evaluation protocol. Different alphabets also define different random-reference distributions.

## Interpretation

- Scores measure compatibility within an assumed helical state; they do not establish folding, binding, membrane insertion, biological activity or experimental stability.
- PDB export builds an ideal helix; it is not structure prediction or energy minimization.
- The 21-point Penetration scan changes an angular exposure mask, with homogeneous endpoints. It is not a free-energy profile along a physical insertion coordinate.
- Calibration and reference samples, orientation and seeds affect results. Finite sampling variation should not be confused with a systematic chemical effect.
- Execution, invariance and objective-consistency tests do not validate the parameter sources or experimental predictions.
