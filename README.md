# Helix Design Studio

[Web server](https://helix.simbioslab.com/) · [Method article](https://doi.org/10.1021/acs.jctc.6c01278) · [Citation](CITATION.cff)

Helix Design Studio scores and designs peptide sequences in a **predefined alpha-helical state**, with polar, apolar and membrane-like environments. Its browser interface and FastAPI API support sequence scoring, single and family design, environmental comparisons, specificity design, exposure scans and idealized PDB export.

**Known issue:** helix-neighbor scores depend on the order of the requested residue alphabet because of a matrix-indexing defect. See [KNOWN_ISSUES.md](KNOWN_ISSUES.md) before interpreting scores or designs.

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

The image runs as an unprivileged user. Set `PORT` to change its internal listening port and adjust the port mapping accordingly. Container execution has not been tested.

## API example and reproducibility

With the server running:

```bash
curl --fail-with-body http://127.0.0.1:8000/api/score \
  -H 'Content-Type: application/json' \
  --data @examples/score.json
```

This example uses small calibration and reference samples for a quick execution check. Scientific applications require sampling sizes appropriate to the question. Record the source commit, complete request, residue alphabet **including its order**, random seeds and dependency versions when reporting results.

## Tests

```bash
python -m pip install -r requirements-dev.txt
python -m unittest discover -s tests -v
```

The tests cover the computational API, streaming responses, static assets, idealized PDB export, invalid inputs and disabled usage logging. The alphabet-order regression test is marked as an expected failure for the defect described in [KNOWN_ISSUES.md](KNOWN_ISSUES.md). They check implementation behavior; they do not establish sampling convergence or experimental validity.

## Configuration and data

| Variable | Default | Purpose |
| --- | --- | --- |
| `HELIX_CACHE_DIR` | `helix-design-studio-cache` in the OS temporary directory | Generated calibration statistics. |
| `HELIX_USAGE_LOGGING` | `0` | Set to `1` to enable local usage-event logging. |
| `HELIX_VAR_DIR` | A directory in the OS temporary directory | Optional usage logs; choose a persistent directory if required. |
| `QFOLD_CODE_DIR` | This repository's `Code/` | Override the location of scoring modules and parameter tables. |
| `PORT` | `8000` | Container listening port. |

With logging disabled, the application creates no usage-event log or persistent browser visitor identifier. When enabled, events include a browser identifier, hashes of IP address, user agent and sequence, request settings, origin/referrer, timing and error information. Hashes do not guarantee anonymity. Events are stored on the server; server/proxy access logs are configured separately.

`GET /api/health` reports service configuration and, when logging is enabled, recent aggregate failures. Calibration caches are disposable; clear them after changing scoring code or parameter tables.

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
