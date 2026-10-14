# API guide

The HTTP API provides the same calculations as the browser interface. Send JSON requests to a server running the version you intend to use, and retain the complete JSON responses for reproducibility.

## Interactive documentation

On a local installation, open [Swagger UI](http://127.0.0.1:8000/docs), [ReDoc](http://127.0.0.1:8000/redoc), or the [OpenAPI schema](http://127.0.0.1:8000/openapi.json). Swagger UI shows request fields and validation limits and lets you execute requests with **Try it out**.

The public server offers [Swagger UI](https://helix.simbioslab.com/docs), [ReDoc](https://helix.simbioslab.com/redoc), and its [OpenAPI schema](https://helix.simbioslab.com/openapi.json). Check the schema's `info.version` and each result's `provenance.version` before comparing calculations from different servers. For a specified release, install its Git tag and run it locally using the [installation instructions](../README.md#installation).

## Available calculations

All calculation endpoints below accept HTTP `POST` with `Content-Type: application/json`.

| Endpoint | Purpose | Main response fields |
|---|---|---|
| `/api/score` | Score one sequence in one environment | `sequence`, `energy`, `z_score`, `breakdown`, `environment`, `calibration`, `random_reference` |
| `/api/design` | Search for one sequence in a target environment | Standard score fields, `design_geometry`, `design_settings` |
| `/api/compare` | Evaluate the same sequence in two environments | `environment_a`, `environment_b`, `energy_gap`, `z_score_gap` |
| `/api/optimal-penetration` | Scan angular membrane exposure from 0 to 100% | `profile`, `best_penetration_percent`, `best_score`, `scan_settings` |
| `/api/design-family` | Generate a sequence family in one environment | `members`, `family_size_returned` |
| `/api/design-family-stream` | Receive family candidates as they become available | Newline-delimited JSON events |
| `/api/cross-design` | Search for compatibility with two environments | `sequence`, `environment_a`, `environment_b`, `cross_design` |
| `/api/specificity-design` | Search for candidates for a target versus off-target environments | `candidates`, `specificity_collection`; top-ranked candidate also exposed as `sequence`, `target_environment`, `off_target_environments`, `specificity_design` |
| `/api/specificity-design-stream` | Receive specificity candidates as they become available | Newline-delimited JSON events |
| `/api/pdb` | Construct an ideal alpha-helical PDB for a sequence | `sequence`, `pdb` |

Calculation responses include `provenance` with the software version, parameter checksum, validated request, seed protocol and runtime versions. Streaming responses attach provenance to the metadata, candidate and final payloads. The streaming content type is `application/x-ndjson`: read one JSON object per line and inspect the `type` field (`meta`, `phase`, `progress`, `candidate`, `final`, or `error`). An error event can occur after the HTTP response has started, so HTTP status alone does not establish that a streamed calculation completed.

For family design, iterate over `members`. For specificity design, iterate over `candidates`, whose entries each contain a sequence, target/off-target scores and `specificity_design` metrics. `specificity_collection.requested` and `.returned` identify the requested and recovered collection sizes; the top-level candidate fields duplicate the highest-ranked entry for convenient inspection. In both streaming workflows, each `candidate` event carries that candidate in `payload.member`, and the `final` event carries the complete collection in `payload`.

## Inputs and scoring conventions

- `sequence` uses standard one-letter amino-acid codes. The `residues` list specifies the allowed alphabet and must include every residue in a supplied sequence. Its order is part of the random sampling protocol.
- `environment.preset` is one of `homogeneous_polar`, `homogeneous_apolar`, `interfacial_neg`, `interfacial_neu`, or `interfacial_pos`. Interfacial presets represent anionic, neutral or cationic membrane surfaces. Their default apolar-sector half-width is 90°. `wheel_halfwidth_deg` can override this angular exposure. `membrane_charge` can override the surface-charge convention with `neg`, `neu`, or `pos`.
- `weights` can override the six Hamiltonian weights. Omitted weights use the release defaults shown in the OpenAPI schema.
- `n_decoys` sets the calibration sample for standardizing Hamiltonian components. `n_random` sets the separate random-sequence reference used by Score and Penetration for the total-score Z-score. Design and Compare use a fixed reference of 500 sequences per environment in this version.
- `seed` makes a request reproducible under the same software, parameters, dependency versions and complete protocol. Calibration, reference sampling and searches use separate random streams.
- `length`, `steps`, and `restarts` control sequence searches. Additional family, cross-environment and specificity options are described in the interactive schema.

The standard interfacial score rotates **each sequence by its own hydrophobic moment**, pointing that moment at the center of the apolar sector. Query sequences, calibration decoys and random reference sequences use the same rule. If the moment is numerically zero, the preset phase is retained. Searches optimize a fixed wheel geometry and additionally report the standard oriented score in their results; `design_geometry` identifies the geometry actually optimized. Lower scores indicate better compatibility with the assumed helical state within a specified protocol. They are dimensionless comparative scores and do not establish antimicrobial activity or membrane selectivity.

The Penetration endpoint scans 21 exposure settings in 5% increments. Zero percent is a homogeneous polar environment; 100% is a homogeneous apolar environment. Intermediate settings use an apolar-sector half-width of `40 + 1.1 × percent` degrees and the same hydrophobic-moment alignment rule. This percentage describes an **angular exposure mask**, not a physical insertion depth.

## Examples with curl

Start the server as described in the [README](../README.md#installation), then run these commands from the repository root. The example files use short searches and only 100 calibration sequences to demonstrate the interface quickly; these settings do not establish sampling convergence for scientific use.

**Score a short sequence:**

```bash
curl --fail-with-body http://127.0.0.1:8000/api/score \
  -H 'Content-Type: application/json' \
  --data @examples/score.json > score-result.json
```

Inspect `energy`, `z_score`, and the six standardized components in `breakdown`. The response's `environment.wheel_phase_deg` records the phase used for this sequence. `calibration.orientation` identifies its orientation protocol.

**Design a 12-residue candidate for an anionic membrane environment:**

```bash
curl --fail-with-body http://127.0.0.1:8000/api/design \
  -H 'Content-Type: application/json' \
  --data @examples/design.json > design-result.json
```

The example allows A, L, K, E, G, F and S and runs one 1,000-step search. Read `sequence` for the candidate, `design_geometry.energy` for the optimized objective, and the top-level `energy` and `z_score` for its standard oriented score.

**Compare magainin 2 in an anionic interface and a polar medium:**

```bash
curl --fail-with-body http://127.0.0.1:8000/api/compare \
  -H 'Content-Type: application/json' \
  --data @examples/compare.json > compare-result.json
```

The response contains the complete scores under `environment_a` and `environment_b`, each with its own calibration and reference distribution. Their absolute differences are reported as `energy_gap` and `z_score_gap`. These are model comparisons; they are not experimental activity or binding measurements.

**Scan angular exposure for the same sequence:**

```bash
curl --fail-with-body http://127.0.0.1:8000/api/optimal-penetration \
  -H 'Content-Type: application/json' \
  --data @examples/penetration.json > penetration-result.json
```

Read the 21 records in `profile` to plot score versus exposure. `best_penetration_percent` minimizes the model Z-score among those settings, and `best_score` contains that setting's complete score. It does not predict a physical insertion depth.

## Python client

The supplied [client](../examples/api_client.py) uses only Python's standard library. It reads the same JSON requests and saves the complete response:

```bash
python examples/api_client.py score --output score-result.json
python examples/api_client.py design --output design-result.json
python examples/api_client.py compare --output compare-result.json
python examples/api_client.py penetration --output penetration-result.json
```

Choose another installation with `--base-url`, for example `--base-url https://helix.simbioslab.com`. Check that installation's software version before using or comparing its numerical results.

## Handling errors

Invalid requests return a non-success HTTP status and a JSON `detail` field identifying the issue. Read the schema for the serving release's limits. A busy server or a calculation exceeding its execution limit also returns an error. The Python client and `curl --fail-with-body` report HTTP failures. For streaming endpoints, inspect every event and require a `final` event before treating the run as complete.
