# Helix Design Studio

[Web server](https://helix.simbioslab.com/) · [Method article](https://doi.org/10.1021/acs.jctc.6c01278) · [Citation](CITATION.cff)

Helix Design Studio scores and designs peptide sequences in a **predefined alpha-helical state**, with polar, apolar and membrane-like environments. Its browser interface and FastAPI API support sequence scoring, single and family design, environmental comparisons, specificity design, exposure scans and idealized PDB export.

## Installation

Use **Python 3.12** and the pinned dependencies. Interactive molecular rendering requires a WebGL-capable browser; the viewer is included in the application. Local execution has been tested on macOS with Python 3.12.

```bash
git clone https://github.com/SIMBIOS-USC/helix-design-studio-web.git
cd helix-design-studio-web
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m uvicorn app.main:app --host 127.0.0.1 --port 8000
```

On Windows, activate the environment with `.venv\Scripts\activate`.

Open **http://127.0.0.1:8000**. The interface includes a Guide and examples. API documentation is at `/docs`; the schema is at `/openapi.json`. The first calculation for a calibration configuration samples random sequences and caches the component statistics for subsequent requests.

## Docker

```bash
docker build -t helix-design-studio .
docker run --rm -p 8000:8000 helix-design-studio
```

The image runs as an unprivileged user and defines a health check. Set `PORT` to change its internal listening port and adjust the port mapping accordingly, `WEB_CONCURRENCY` to run several worker processes, and `FORWARDED_ALLOW_IPS` to the address of the reverse proxy whose `X-Forwarded-For` header should be trusted. Calibration caches and optional usage logs are written to `/tmp/helix-cache` and `/tmp/helix-usage`; mount volumes there (for example `-v helix-cache:/tmp/helix-cache`) to keep them across container restarts. Container execution has not been tested.

## API example and reproducibility

With the server running:

```bash
curl --fail-with-body http://127.0.0.1:8000/api/score \
  -H 'Content-Type: application/json' \
  --data @examples/score.json
```

This example uses small calibration and reference samples for a quick execution check. Scientific applications require sampling sizes appropriate to the question. Record the source commit, complete request, residue alphabet (its order determines the random streams used for calibration and search), random seeds and dependency versions when reporting results.

## Model scope and limitations

- Scores measure compatibility with an assumed alpha-helical state. They do not establish folding, binding, membrane insertion, biological activity or experimental stability.
- In interfacial environments every evaluated sequence (query, calibration decoys, random references and search proposals) is oriented with its own hydrophobic moment, so all of them are compared under the same rule. The effective wheel phase is reported with each result.
- Helix-pair propensities *P* enter the energy as −ln *P*; i,i+3 and i,i+4 tables are energies (negative values are stabilizing). All tables are indexed by residue label.
- The membrane-charge term acts only on water-facing residues of charged interfaces; it is zero in neutral interfaces and homogeneous media.
- PDB export constructs an ideal helix; it does not predict or minimize a structure.
- The Penetration scan changes an angular exposure mask. Its coordinate is not a physical insertion depth or a free-energy reaction coordinate.
- Calibration and reference distributions depend on sequence length, alphabet, environment, sample sizes and seeds. Scores from different protocols are not directly interchangeable.

## Tests

```bash
python -m pip install -r requirements-dev.txt
python -m unittest discover -s tests -v
```

The tests cover the computational API, streaming responses, static assets, idealized PDB export, input validation, request limits, usage logging, calibration caching, label-based parameter indexing, literature residue scales and consistency between workflows. They check implementation behavior; they do not establish sampling convergence or experimental validity.

## Configuration and data

| Variable | Default | Purpose |
| --- | --- | --- |
| `HELIX_CACHE_DIR` | `helix-design-studio-cache` in the OS temporary directory | Generated calibration statistics (private to the service user). |
| `HELIX_CACHE_MAX_FILES` | `2000` | Maximum number of cached calibrations; the least recently used are removed. |
| `HELIX_MAX_CONCURRENT_JOBS` | `4` | Calculations run at the same time per server process; further requests receive HTTP 503. |
| `HELIX_MAX_JOB_SECONDS` | `1800` | Time limit for a single calculation. |
| `HELIX_CORS_ORIGINS` | `https://simbios.usc.es` | Comma-separated origins allowed to call the API from other sites; empty disables CORS. |
| `HELIX_USAGE_LOGGING` | `0` | Set to `1` to enable local usage-event logging. |
| `HELIX_VAR_DIR` | A directory in the OS temporary directory | Optional usage logs; choose a persistent directory if required. |
| `HELIX_LOG_HASH_KEY` | Random per process | Secret key for the hashes stored in usage logs; set it to keep hashes comparable across restarts. |
| `HELIX_CODE_DIR` | This repository's `Code/` | Override the location of scoring modules and parameter tables (`QFOLD_CODE_DIR` is accepted as a legacy name). |
| `PORT` | `8000` | Container listening port. |

With logging disabled, the application creates no usage-event log or persistent browser visitor identifier. The visit and run counters shown when logging is enabled are informational: they are reported by browsers and are not protected against inflation. When enabled, events include a random browser identifier, keyed hashes (HMAC-SHA256) of IP address, user agent and sequences, request settings, truncated origin/referrer, timing, HTTP status and error messages. Hashes do not guarantee anonymity. Events are stored on the server; server/proxy access logs are configured separately.

`GET /api/health` reports service configuration and, when logging is enabled, recent server-side failures (HTTP 5xx); it returns 503 while such failures are present in the requested window. Calibration caches are keyed by the parameter tables and settings; they are disposable and can be cleared at any time.

For a public deployment, run the service behind a reverse proxy that terminates TLS and applies per-client rate limits, start uvicorn with `--proxy-headers --forwarded-allow-ips=<proxy address>` so that client addresses are taken from the trusted proxy, and use `--workers` to serve several calculations in parallel.

## Source layout

- `app/`: API, search workflows and static browser interface.
- `Code/`: scoring, calibration and four parameter tables.
- `tests/` and `examples/`: automated checks and an API request example.
- `LICENSES/`, `THIRD_PARTY_NOTICES.md`, `CITATION.cff`: licenses, attribution and citation metadata.

## Citation and license

Please cite the **software version or commit** used and the underlying method:

Daniel Conde-Torres, Rebeca García-Fandiño and Ángel Piñeiro. *Environment-Conditioned Design of α-Helical Peptides*. Journal of Chemical Theory and Computation **22** (18), 9776–9789 (2026). [doi:10.1021/acs.jctc.6c01278](https://doi.org/10.1021/acs.jctc.6c01278).

Software contributors are listed in [CITATION.cff](CITATION.cff). Original project code is distributed under the [MIT license](LICENSE). Included components retain their own terms; see [third-party notices](THIRD_PARTY_NOTICES.md). Institutional names and logos are excluded from the MIT grant. Citation is requested as scholarly practice, not as an additional license condition.

Report problems through the [issue tracker](https://github.com/SIMBIOS-USC/helix-design-studio-web/issues), including the commit, request parameters and steps to reproduce. Omit private sequences and usage logs from public reports.
