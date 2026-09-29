# Helix Design Studio

[Web server](https://helix.simbioslab.com/) · [Method article](https://doi.org/10.1021/acs.jctc.6c01278) · [Citation](CITATION.cff)

Helix Design Studio scores and designs peptide sequences in a **predefined alpha-helical state**, with polar, apolar and membrane-like environments. It provides sequence scoring, single and family design, comparison between environments, specificity design, penetration scans and idealized PDB export through a browser and a FastAPI API.

**Scientific status (0.2.0rc1):** this review candidate corrects residue-label mapping and preserves search/display calibration. Two malformed parameter-table headers have been reconstructed provisionally; their source/order still needs author confirmation. See [PARAMETER_PROVENANCE.md](PARAMETER_PROVENANCE.md) and [KNOWN_ISSUES.md](KNOWN_ISSUES.md). Passing software tests establishes implementation consistency, not scientific validity. This candidate is not yet a validated submission release.

## Run locally

Python **3.12** is the tested interpreter. No database, GPU, quantum software, Node.js or external graphics program is required. The molecular viewer is included locally.

```bash
git clone https://github.com/SIMBIOS-USC/helix-design-studio-web.git
cd helix-design-studio-web
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python -m uvicorn app.main:app --host 127.0.0.1 --port 8000
```

Open **http://127.0.0.1:8000**. The interface includes a Guide and examples. API documentation is at `/docs`; the schema is at `/openapi.json`. The first calculation for a new calibration configuration can take longer because it samples reference sequences and caches their statistics.

On Windows, activate the environment with `.venv\Scripts\activate`. Cross-platform installation and container instructions are provided; the initial local verification was performed on macOS with Python 3.12.

## Run in Docker

```bash
docker build -t helix-design-studio .
docker run --rm -p 8000:8000 helix-design-studio
```

The image runs as an unprivileged user and serves port 8000 by default. `PORT` changes its internal listening port; adjust the port mapping accordingly. Container execution must be verified on the deployment host.

## Reproduce a small example

With the server running:

```bash
curl --fail-with-body http://127.0.0.1:8000/api/score \
  -H 'Content-Type: application/json' \
  --data @examples/score.json
```

The example uses small calibration and reference samples for a quick execution check. It is not a benchmark or a recommended scientific sampling protocol. Record the source commit, request parameters, residue alphabet **including its order**, seed and dependency versions when reporting results. Responses include calibration/reference seeds and an orientation-mode field. Score/Compare automatically align interfacial sequences to the hydrophobic moment; search outputs retain fixed search geometry. A subsequent Score request is a distinct evaluation protocol.

## Verification

```bash
python -m pip install -r requirements-dev.txt
python -m unittest discover -s tests -v
```

The tests exercise the computational API, streaming responses, static assets, idealized PDB export, invalid-input handling and disabled usage logging. Regression tests check raw-score invariance to alphabet order, matrix validation and agreement between search objectives and displayed scores. Small test inputs check execution and consistency, not scientific convergence or performance at the web limits.

## Configuration and data

| Variable | Default | Purpose |
| --- | --- | --- |
| `HELIX_CACHE_DIR` | `helix-design-studio-cache` in the OS temporary directory | Generated calibration statistics. |
| `HELIX_USAGE_LOGGING` | `0` | Set to `1` to enable local usage-event logging. |
| `HELIX_VAR_DIR` | A directory in the OS temporary directory | Location of optional usage logs; choose a persistent directory if needed. |
| `QFOLD_CODE_DIR` | This repository's `Code/` | Optional explicit runtime override; normally unnecessary. |
| `PORT` | `8000` | Container listening port. |

With logging disabled, no application usage-event log or persistent browser visitor identifier is created. When enabled, events include a browser identifier, hashes of IP address, user agent and sequence, request settings, origin/referrer, timing and error information. Hashes are not a guarantee of anonymity. These events remain local to the server; no external analytics service is used. Server/proxy access logs are configured separately.

`GET /api/health` reports service configuration and, when logging is enabled, recent aggregate failures. It is a service check, not a scientific validation check. Calibration caches are disposable; clear them after changes to scoring code or tables. No user logs, precomputed caches or generated results are distributed here.

## Source layout

- `app/`: API, classical search workflows and static browser interface.
- `Code/`: scoring/calibration routines and the four runtime parameter tables.
- `tests/` and `examples/`: small reproducible checks and an API example.
- `LICENSES/`, `THIRD_PARTY_NOTICES.md`, `CITATION.cff`: licenses, attribution and citation metadata.

This distribution excludes historical quantum backends, table-generation experiments, research datasets, manuscripts, private deployment scripts and source repository history. The production website is maintained separately; its running commit has not been verified against this snapshot. Changes introduced for this distribution are described in [PROVENANCE.md](PROVENANCE.md).

## Citation and reuse

If this software contributes to your research, please cite **the exact software version/commit** and the underlying method:

Daniel Conde-Torres, Rebeca García-Fandiño and Ángel Piñeiro. *Environment-Conditioned Design of α-Helical Peptides*. Journal of Chemical Theory and Computation **22** (18), 9776–9789 (2026). [doi:10.1021/acs.jctc.6c01278](https://doi.org/10.1021/acs.jctc.6c01278).

Software contributors are listed in [CITATION.cff](CITATION.cff). A software-paper citation and an archived software-release DOI can be added when available; neither is assigned by this repository. Citation is requested as scholarly practice, not as an additional restriction on the MIT license.

Original project code is distributed under the [MIT license](LICENSE). Components carrying other licenses retain their own terms; see [third-party notices](THIRD_PARTY_NOTICES.md). Institutional names and logos identify the project and are excluded from the MIT grant.

Report reproducible problems through this repository's issue tracker, including the commit and request parameters. Do not include private sequences or usage logs in public reports.
