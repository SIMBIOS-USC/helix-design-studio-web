from __future__ import annotations

import hashlib
import json
import os
import tempfile
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

from fastapi import FastAPI, HTTPException, Request
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse, StreamingResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel, Field, field_validator

from .engine import (
    DEFAULT_ALPHABET_16,
    DEFAULT_WEIGHTS,
    ENVIRONMENT_PRESETS,
    compare_environments,
    cross_environment_design,
    design_specificity,
    design_specificity_stream,
    design_family,
    design_family_stream,
    design_sequence,
    generate_helix_pdb,
    optimize_penetration,
    score_sequence,
)

APP_DIR = Path(__file__).resolve().parent
STATIC_DIR = APP_DIR / "static"
USAGE_LOGGING_ENABLED = os.environ.get("HELIX_USAGE_LOGGING", "0") == "1"
VAR_DIR = Path(
    os.environ.get("HELIX_VAR_DIR")
    or Path(tempfile.gettempdir()) / "helix-design-studio-usage"
).expanduser()
USAGE_LOG_PATH = VAR_DIR / "usage_events.jsonl"

app = FastAPI(title="Helix Design Studio", version="0.2.0rc1")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["https://simbios.usc.es"],
    allow_methods=["GET", "POST", "OPTIONS"],
    allow_headers=["*"],
)

app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")

MAX_SEQUENCE_LENGTH = 40
MAX_ALLOWED_RESIDUES = 20
MIN_ALLOWED_RESIDUES = 2
MAX_RANDOM_REFERENCE = 2000
MAX_N_DECOYS = 10000
MAX_DESIGN_STEPS = 50000
MAX_FAMILY_STEPS = 20000
MAX_CROSS_STEPS = 60000
MAX_RESTARTS = 20
MAX_FAMILY_SIZE = 12
MAX_SPECIFICITY_CANDIDATES = 10
MAX_OVERSAMPLE_FACTOR = 5
MAX_PENETRATION_SCAN_POINTS = 21
SIMBIOS_CONTACT_MSG = "For heavier runs or custom large-scale calculations, please contact SIMBIOS via https://simbios.usc.es/."
FULL_ALPHABET_20 = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]


def _clean_sequence(value: str) -> str:
    seq = str(value or "").strip().upper()
    if not seq:
        raise ValueError("Sequence cannot be empty.")
    if len(seq) > MAX_SEQUENCE_LENGTH:
        raise ValueError(
            f"Sequence length exceeds the current web limit of {MAX_SEQUENCE_LENGTH} residues. {SIMBIOS_CONTACT_MSG}"
        )
    return seq


def _validate_residues(value: List[str]) -> List[str]:
    residues = [str(aa).strip().upper() for aa in value if str(aa).strip()]
    unique: List[str] = []
    for aa in residues:
        if aa not in unique:
            unique.append(aa)
    if len(unique) < MIN_ALLOWED_RESIDUES:
        raise ValueError("Please provide at least two distinct allowed residue types.")
    if len(unique) > MAX_ALLOWED_RESIDUES:
        raise ValueError(
            f"Allowed-residue sets are limited to {MAX_ALLOWED_RESIDUES} residue types in the web app. {SIMBIOS_CONTACT_MSG}"
        )
    return unique


def _hash_value(value: str) -> str:
    return hashlib.sha256(value.encode("utf-8")).hexdigest()[:16]


def _sequence_summary(sequence: Optional[str]) -> Dict[str, Any]:
    if not sequence:
        return {}
    seq = str(sequence).strip().upper()
    if not seq:
        return {}
    return {
        "sequence_length": len(seq),
        "sequence_hash": _hash_value(seq),
    }


def _environment_summary(environment: Optional[Dict[str, Any]], prefix: str = "environment") -> Dict[str, Any]:
    if not environment:
        return {}
    return {
        f"{prefix}_preset": environment.get("preset"),
        f"{prefix}_membrane_charge": environment.get("membrane_charge"),
        f"{prefix}_wheel_halfwidth_deg": environment.get("wheel_halfwidth_deg"),
    }


def _client_ip_from_request(request: Request) -> str:
    forwarded = request.headers.get("x-forwarded-for")
    if forwarded:
        return forwarded.split(",")[0].strip()
    if request.client and request.client.host:
        return request.client.host
    return "unknown"


def _usage_context(request: Request) -> Dict[str, Any]:
    anonymous_id = request.headers.get("x-anonymous-id", "").strip() or None
    user_agent = request.headers.get("user-agent", "").strip() or "unknown"
    client_ip = _client_ip_from_request(request)
    return {
        "anonymous_user_id": anonymous_id,
        "user_agent_hash": _hash_value(user_agent),
        "client_ip_hash": _hash_value(client_ip),
        "origin": request.headers.get("origin"),
        "referer": request.headers.get("referer"),
    }


def _request_summary(payload: Any) -> Dict[str, Any]:
    data = payload.model_dump() if hasattr(payload, "model_dump") else dict(payload)
    summary: Dict[str, Any] = {}
    summary.update(_sequence_summary(data.get("sequence")))
    if "length" in data:
        summary["length"] = data["length"]
    if "residues" in data:
        summary["alphabet_size"] = len(data["residues"])
    if "steps" in data:
        summary["steps"] = data["steps"]
    if "restarts" in data:
        summary["restarts"] = data["restarts"]
    if "family_size" in data:
        summary["family_size"] = data["family_size"]
    if "num_sequences" in data:
        summary["num_sequences"] = data["num_sequences"]
    if "n_random" in data:
        summary["n_random"] = data["n_random"]
    if "lambda_gap" in data:
        summary["lambda_gap"] = data["lambda_gap"]
    if "lambda_balance" in data:
        summary["lambda_balance"] = data["lambda_balance"]
    if "percent_min" in data:
        summary["percent_min"] = data["percent_min"]
    if "percent_max" in data:
        summary["percent_max"] = data["percent_max"]
    if "percent_step" in data:
        summary["percent_step"] = data["percent_step"]
    if "environment" in data:
        summary.update(_environment_summary(data.get("environment"), "environment"))
    if "environment_a" in data:
        summary.update(_environment_summary(data.get("environment_a"), "environment_a"))
    if "environment_b" in data:
        summary.update(_environment_summary(data.get("environment_b"), "environment_b"))
    return summary


def _append_usage_event(event: Dict[str, Any]) -> None:
    if not USAGE_LOGGING_ENABLED:
        return
    VAR_DIR.mkdir(mode=0o700, parents=True, exist_ok=True)
    descriptor = os.open(USAGE_LOG_PATH, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o600)
    with os.fdopen(descriptor, "a", encoding="utf-8") as handle:
        handle.write(json.dumps(event, ensure_ascii=True) + "\n")


def _count_usage_events(endpoint: Optional[str] = None, status: Optional[str] = None) -> int:
    if not USAGE_LOGGING_ENABLED or not USAGE_LOG_PATH.exists():
        return 0
    count = 0
    with USAGE_LOG_PATH.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            if endpoint is not None and event.get("endpoint") != endpoint:
                continue
            if status is not None and event.get("status") != status:
                continue
            count += 1
    return count


def _usage_metrics() -> Dict[str, int]:
    if not USAGE_LOGGING_ENABLED or not USAGE_LOG_PATH.exists():
        return {"total_visits": 0, "unique_visitors": 0, "total_runs": 0}

    total_visits = 0
    unique_visitors: set[str] = set()
    total_runs = 0

    with USAGE_LOG_PATH.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue
            if event.get("status") != "ok":
                continue

            endpoint = event.get("endpoint")
            anonymous_id = event.get("anonymous_user_id")
            if endpoint == "visit":
                total_visits += 1
                if anonymous_id:
                    unique_visitors.add(str(anonymous_id))
                continue

            total_runs += 1

    return {
        "total_visits": total_visits,
        "unique_visitors": len(unique_visitors),
        "total_runs": total_runs,
    }


DEFAULT_HEALTH_WINDOW_MINUTES = 60


def _recent_usage_summary(window_minutes: int) -> Dict[str, Any]:
    """Resume los eventos de uso de la ventana reciente, agrupados por endpoint.

    Se usa en /api/health: los endpoints de streaming capturan sus excepciones y
    las emiten dentro del stream, por lo que HTTP siempre responde 200 y un
    monitor externo basado en codigos de estado no puede detectar los fallos.
    """
    summary: Dict[str, Any] = {
        "window_minutes": window_minutes,
        "events": 0,
        "errors": 0,
        "failing_endpoints": {},
        "last_error": None,
    }
    if not USAGE_LOGGING_ENABLED or not USAGE_LOG_PATH.exists():
        return summary

    cutoff = datetime.now(timezone.utc) - timedelta(minutes=window_minutes)
    failing: Dict[str, int] = {}
    last_error: Optional[Dict[str, Any]] = None

    try:
        handle = USAGE_LOG_PATH.open("r", encoding="utf-8")
    except OSError:
        summary["read_error"] = "Usage log is unavailable."
        return summary

    with handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                event = json.loads(line)
            except json.JSONDecodeError:
                continue

            raw_timestamp = event.get("timestamp_utc")
            if not raw_timestamp:
                continue
            try:
                stamp = datetime.fromisoformat(str(raw_timestamp))
            except ValueError:
                continue
            if stamp.tzinfo is None:
                stamp = stamp.replace(tzinfo=timezone.utc)
            if stamp < cutoff:
                continue

            summary["events"] += 1
            if event.get("status") != "error":
                continue

            summary["errors"] += 1
            endpoint = str(event.get("endpoint", "unknown"))
            failing[endpoint] = failing.get(endpoint, 0) + 1
            last_error = {
                "timestamp_utc": raw_timestamp,
                "endpoint": endpoint,
            }

    summary["failing_endpoints"] = dict(sorted(failing.items(), key=lambda kv: (-kv[1], kv[0])))
    summary["last_error"] = last_error
    return summary


def _log_usage_event(
    request: Request,
    endpoint: str,
    payload: Any,
    status: str,
    duration_s: float,
    extra: Optional[Dict[str, Any]] = None,
) -> None:
    if not USAGE_LOGGING_ENABLED:
        return
    event = {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "endpoint": endpoint,
        "status": status,
        "duration_s": round(duration_s, 4),
        **_usage_context(request),
        **_request_summary(payload),
    }
    if extra:
        event.update(extra)
    _append_usage_event(event)


class EnvironmentConfig(BaseModel):
    preset: str = "interfacial_neg"
    wheel_phase_deg: Optional[float] = None
    wheel_halfwidth_deg: Optional[float] = None
    membrane_charge: Optional[str] = None


class WeightsConfig(BaseModel):
    lambda_env: float = DEFAULT_WEIGHTS["lambda_env"]
    lambda_charge: float = DEFAULT_WEIGHTS["lambda_charge"]
    lambda_local: float = DEFAULT_WEIGHTS["lambda_local"]
    lambda_pairwise: float = DEFAULT_WEIGHTS["lambda_pairwise"]
    lambda_helix_pairs: float = DEFAULT_WEIGHTS["lambda_helix_pairs"]
    lambda_electrostatic: float = DEFAULT_WEIGHTS["lambda_electrostatic"]


class ScoreRequest(BaseModel):
    sequence: str
    residues: List[str] = Field(default_factory=lambda: list(FULL_ALPHABET_20))
    environment: EnvironmentConfig = Field(default_factory=EnvironmentConfig)
    weights: WeightsConfig = Field(default_factory=WeightsConfig)
    n_decoys: int = Field(default=2000, ge=100, le=MAX_N_DECOYS)
    n_random: int = Field(default=500, ge=100, le=MAX_RANDOM_REFERENCE)
    seed: int = 42

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, value: str) -> str:
        return _clean_sequence(value)

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, value: List[str]) -> List[str]:
        return _validate_residues(value)


class DesignRequest(BaseModel):
    length: int = Field(default=20, ge=4, le=MAX_SEQUENCE_LENGTH)
    residues: List[str] = Field(default_factory=lambda: list(DEFAULT_ALPHABET_16))
    environment: EnvironmentConfig = Field(default_factory=EnvironmentConfig)
    weights: WeightsConfig = Field(default_factory=WeightsConfig)
    n_decoys: int = Field(default=2000, ge=100, le=MAX_N_DECOYS)
    steps: int = Field(default=12000, ge=1000, le=MAX_DESIGN_STEPS)
    restarts: int = Field(default=8, ge=1, le=MAX_RESTARTS)
    seed: int = 42

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, value: List[str]) -> List[str]:
        return _validate_residues(value)


class CompareRequest(BaseModel):
    sequence: str
    residues: List[str] = Field(default_factory=lambda: list(FULL_ALPHABET_20))
    environment_a: EnvironmentConfig
    environment_b: EnvironmentConfig
    weights: WeightsConfig = Field(default_factory=WeightsConfig)
    n_decoys: int = Field(default=2000, ge=100, le=MAX_N_DECOYS)
    seed: int = 42

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, value: str) -> str:
        return _clean_sequence(value)

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, value: List[str]) -> List[str]:
        return _validate_residues(value)


class FamilyDesignRequest(BaseModel):
    length: int = Field(default=20, ge=4, le=MAX_SEQUENCE_LENGTH)
    residues: List[str] = Field(default_factory=lambda: list(DEFAULT_ALPHABET_16))
    environment: EnvironmentConfig = Field(default_factory=EnvironmentConfig)
    weights: WeightsConfig = Field(default_factory=WeightsConfig)
    n_decoys: int = Field(default=2000, ge=100, le=MAX_N_DECOYS)
    steps: int = Field(default=10000, ge=1000, le=MAX_FAMILY_STEPS)
    restarts: int = Field(default=6, ge=1, le=MAX_RESTARTS)
    family_size: int = Field(default=8, ge=2, le=MAX_FAMILY_SIZE)
    oversample_factor: int = Field(default=3, ge=1, le=MAX_OVERSAMPLE_FACTOR)
    seed: int = 42

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, value: List[str]) -> List[str]:
        return _validate_residues(value)


class CrossDesignRequest(BaseModel):
    length: int = Field(default=20, ge=4, le=MAX_SEQUENCE_LENGTH)
    residues: List[str] = Field(default_factory=lambda: list(DEFAULT_ALPHABET_16))
    environment_a: EnvironmentConfig
    environment_b: EnvironmentConfig
    lambda_gap: float = Field(default=0.5, ge=0.0, le=10.0)
    weights: WeightsConfig = Field(default_factory=WeightsConfig)
    n_decoys: int = Field(default=2000, ge=100, le=MAX_N_DECOYS)
    steps: int = Field(default=16000, ge=1000, le=MAX_CROSS_STEPS)
    restarts: int = Field(default=10, ge=1, le=MAX_RESTARTS)
    seed: int = 42

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, value: List[str]) -> List[str]:
        return _validate_residues(value)


class SpecificityDesignRequest(BaseModel):
    length: int = Field(default=20, ge=4, le=MAX_SEQUENCE_LENGTH)
    residues: List[str] = Field(default_factory=lambda: list(FULL_ALPHABET_20))
    target_environment: EnvironmentConfig
    off_target_environments: List[EnvironmentConfig]
    num_sequences: int = Field(default=10, ge=1, le=MAX_SPECIFICITY_CANDIDATES)
    lambda_balance: float = Field(default=0.0, ge=0.0, le=1.0)
    weights: WeightsConfig = Field(default_factory=WeightsConfig)
    n_decoys: int = Field(default=10000, ge=100, le=MAX_N_DECOYS)
    steps: int = Field(default=60000, ge=1000, le=MAX_CROSS_STEPS)
    restarts: int = Field(default=20, ge=1, le=MAX_RESTARTS)
    seed: int = 42

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, value: List[str]) -> List[str]:
        return _validate_residues(value)

    @field_validator("off_target_environments")
    @classmethod
    def validate_off_targets(cls, value: List[EnvironmentConfig]) -> List[EnvironmentConfig]:
        if not value:
            raise ValueError("Please provide at least one off-target environment.")
        return value

    @field_validator("seed")
    @classmethod
    def validate_target_not_in_off_targets(cls, value: int, info: Any) -> int:
        data = getattr(info, "data", {}) or {}
        target = data.get("target_environment")
        off_targets = data.get("off_target_environments") or []
        if not target:
            return value

        def env_key(env: Any) -> tuple[Any, Any, Any, Any]:
            if hasattr(env, "preset"):
                return (
                    getattr(env, "preset", None),
                    getattr(env, "wheel_phase_deg", None),
                    getattr(env, "wheel_halfwidth_deg", None),
                    getattr(env, "membrane_charge", None),
                )
            return (
                env.get("preset"),
                env.get("wheel_phase_deg"),
                env.get("wheel_halfwidth_deg"),
                env.get("membrane_charge"),
            )

        target_key = env_key(target)
        if any(env_key(env) == target_key for env in off_targets):
            raise ValueError("The target environment cannot also be listed as an off-target.")
        return value


class PenetrationOptimizationRequest(BaseModel):
    sequence: str
    residues: List[str] = Field(default_factory=lambda: list(FULL_ALPHABET_20))
    environment: EnvironmentConfig
    weights: WeightsConfig = Field(default_factory=WeightsConfig)
    n_decoys: int = Field(default=2000, ge=100, le=MAX_N_DECOYS)
    n_random: int = Field(default=500, ge=100, le=MAX_RANDOM_REFERENCE)
    percent_min: int = Field(default=0, ge=0, le=100)
    percent_max: int = Field(default=100, ge=0, le=100)
    percent_step: int = Field(default=5, ge=1, le=20)
    seed: int = 42

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, value: str) -> str:
        return _clean_sequence(value)

    @field_validator("residues")
    @classmethod
    def validate_residues(cls, value: List[str]) -> List[str]:
        return _validate_residues(value)

    @field_validator("environment")
    @classmethod
    def validate_environment(cls, value: EnvironmentConfig) -> EnvironmentConfig:
        if not str(value.preset).startswith("interfacial_"):
            raise ValueError("Optimal penetration can only be scanned for interfacial presets.")
        return value

    @field_validator("seed")
    @classmethod
    def validate_scan_geometry(cls, value: int, info: Any) -> int:
        data = getattr(info, "data", {}) or {}
        percent_min = int(data.get("percent_min", 0))
        percent_max = int(data.get("percent_max", 100))
        percent_step = int(data.get("percent_step", 5))
        if percent_min != 0 or percent_max != 100 or percent_step != 5:
            raise ValueError("Web penetration scans are fixed to 0-100% in 5% increments.")
        n_points = ((percent_max - percent_min) // percent_step) + 1
        if n_points > MAX_PENETRATION_SCAN_POINTS:
            raise ValueError(
                f"Penetration scans are limited to {MAX_PENETRATION_SCAN_POINTS} sampled depths in the web app. {SIMBIOS_CONTACT_MSG}"
            )
        return value


class PdbRequest(BaseModel):
    sequence: str

    @field_validator("sequence")
    @classmethod
    def validate_sequence(cls, value: str) -> str:
        return _clean_sequence(value)


@app.get("/")
def index() -> FileResponse:
    return FileResponse(STATIC_DIR / "index.html")


@app.get("/api/health")
def health(window_minutes: int = DEFAULT_HEALTH_WINDOW_MINUTES) -> JSONResponse:
    window_minutes = max(1, min(int(window_minutes), 60 * 24 * 7))
    recent = _recent_usage_summary(window_minutes)
    degraded = recent["errors"] > 0

    body: Dict[str, Any] = {
        "status": "degraded" if degraded else "ok",
        "service": "Helix Design Studio",
        "usage_logging": USAGE_LOGGING_ENABLED,
        "recent": recent,
        "alphabet_default": DEFAULT_ALPHABET_16,
        "environment_presets": ENVIRONMENT_PRESETS,
        "default_weights": DEFAULT_WEIGHTS,
    }
    # 503 para que un monitor externo basado en codigos HTTP vea el fallo: los
    # endpoints de streaming responden 200 aunque la ejecucion falle.
    return JSONResponse(content=body, status_code=503 if degraded else 200)


@app.post("/api/visit")
def register_visit(request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        _log_usage_event(request, "visit", {}, "ok", time.perf_counter() - started)
        return _usage_metrics()
    except Exception as exc:
        _log_usage_event(request, "visit", {}, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/api/usage-metrics")
def usage_metrics() -> Dict[str, Any]:
    return _usage_metrics()


@app.post("/api/score")
def score(payload: ScoreRequest, request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        result = score_sequence(
            sequence=payload.sequence,
            residues=payload.residues,
            environment=payload.environment.model_dump(),
            weights=payload.weights.model_dump(),
            n_decoys=payload.n_decoys,
            n_random=payload.n_random,
            seed=payload.seed,
        )
        _log_usage_event(request, "score", payload, "ok", time.perf_counter() - started)
        return result
    except Exception as exc:
        _log_usage_event(request, "score", payload, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/design")
def design(payload: DesignRequest, request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        result = design_sequence(
            length=payload.length,
            residues=payload.residues,
            environment=payload.environment.model_dump(),
            weights=payload.weights.model_dump(),
            n_decoys=payload.n_decoys,
            steps=payload.steps,
            restarts=payload.restarts,
            seed=payload.seed,
        )
        _log_usage_event(request, "design", payload, "ok", time.perf_counter() - started)
        return result
    except Exception as exc:
        _log_usage_event(request, "design", payload, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/compare")
def compare(payload: CompareRequest, request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        result = compare_environments(
            sequence=payload.sequence,
            residues=payload.residues,
            environment_a=payload.environment_a.model_dump(),
            environment_b=payload.environment_b.model_dump(),
            weights=payload.weights.model_dump(),
            n_decoys=payload.n_decoys,
            seed=payload.seed,
        )
        _log_usage_event(request, "compare", payload, "ok", time.perf_counter() - started)
        return result
    except Exception as exc:
        _log_usage_event(request, "compare", payload, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/design-family")
def family_design(payload: FamilyDesignRequest, request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        result = design_family(
            length=payload.length,
            residues=payload.residues,
            environment=payload.environment.model_dump(),
            weights=payload.weights.model_dump(),
            n_decoys=payload.n_decoys,
            steps=payload.steps,
            restarts=payload.restarts,
            family_size=payload.family_size,
            oversample_factor=payload.oversample_factor,
            seed=payload.seed,
        )
        _log_usage_event(
            request,
            "design-family",
            payload,
            "ok",
            time.perf_counter() - started,
            {"family_size_returned": result.get("family_size_returned")},
        )
        return result
    except Exception as exc:
        _log_usage_event(request, "design-family", payload, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/design-family-stream")
def family_design_stream_route(payload: FamilyDesignRequest, request: Request) -> StreamingResponse:
    def generate() -> Any:
        started = time.perf_counter()
        try:
            final_payload: Optional[Dict[str, Any]] = None
            for update in design_family_stream(
                length=payload.length,
                residues=payload.residues,
                environment=payload.environment.model_dump(),
                weights=payload.weights.model_dump(),
                n_decoys=payload.n_decoys,
                steps=payload.steps,
                restarts=payload.restarts,
                family_size=payload.family_size,
                oversample_factor=payload.oversample_factor,
                seed=payload.seed,
            ):
                if update.get("type") == "final":
                    final_payload = update.get("payload", {})
                yield json.dumps(update) + "\n"
            _log_usage_event(
                request,
                "design-family-stream",
                payload,
                "ok",
                time.perf_counter() - started,
                {"family_size_returned": (final_payload or {}).get("family_size_returned")},
            )
        except Exception as exc:
            _log_usage_event(
                request,
                "design-family-stream",
                payload,
                "error",
                time.perf_counter() - started,
                {"error": str(exc)},
            )
            yield json.dumps({"type": "error", "detail": str(exc)}) + "\n"

    return StreamingResponse(generate(), media_type="application/x-ndjson")


@app.post("/api/cross-design")
def cross_design(payload: CrossDesignRequest, request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        result = cross_environment_design(
            length=payload.length,
            residues=payload.residues,
            environment_a=payload.environment_a.model_dump(),
            environment_b=payload.environment_b.model_dump(),
            lambda_gap=payload.lambda_gap,
            weights=payload.weights.model_dump(),
            n_decoys=payload.n_decoys,
            steps=payload.steps,
            restarts=payload.restarts,
            seed=payload.seed,
        )
        _log_usage_event(request, "cross-design", payload, "ok", time.perf_counter() - started)
        return result
    except Exception as exc:
        _log_usage_event(request, "cross-design", payload, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/specificity-design")
def specificity_design(payload: SpecificityDesignRequest, request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        result = design_specificity(
            length=payload.length,
            residues=payload.residues,
            target_environment=payload.target_environment.model_dump(),
            off_target_environments=[env.model_dump() for env in payload.off_target_environments],
            num_sequences=payload.num_sequences,
            lambda_balance=payload.lambda_balance,
            weights=payload.weights.model_dump(),
            n_decoys=payload.n_decoys,
            steps=payload.steps,
            restarts=payload.restarts,
            seed=payload.seed,
        )
        _log_usage_event(request, "specificity-design", payload, "ok", time.perf_counter() - started)
        return result
    except Exception as exc:
        _log_usage_event(request, "specificity-design", payload, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/specificity-design-stream")
def specificity_design_stream_route(payload: SpecificityDesignRequest, request: Request) -> StreamingResponse:
    def generate() -> Any:
        started = time.perf_counter()
        try:
            final_payload: Optional[Dict[str, Any]] = None
            for update in design_specificity_stream(
                length=payload.length,
                residues=payload.residues,
                target_environment=payload.target_environment.model_dump(),
                off_target_environments=[env.model_dump() for env in payload.off_target_environments],
                num_sequences=payload.num_sequences,
                lambda_balance=payload.lambda_balance,
                weights=payload.weights.model_dump(),
                n_decoys=payload.n_decoys,
                steps=payload.steps,
                restarts=payload.restarts,
                seed=payload.seed,
            ):
                if update.get("type") == "final":
                    final_payload = update.get("payload", {})
                yield json.dumps(update) + "\n"
            _log_usage_event(
                request,
                "specificity-design-stream",
                payload,
                "ok",
                time.perf_counter() - started,
                {"final_sequence": (final_payload or {}).get("sequence")},
            )
        except Exception as exc:
            _log_usage_event(
                request,
                "specificity-design-stream",
                payload,
                "error",
                time.perf_counter() - started,
                {"error": str(exc)},
            )
            yield json.dumps({"type": "error", "detail": str(exc)}) + "\n"

    return StreamingResponse(generate(), media_type="application/x-ndjson")


@app.post("/api/optimal-penetration")
def optimal_penetration(payload: PenetrationOptimizationRequest, request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        result = optimize_penetration(
            sequence=payload.sequence,
            residues=payload.residues,
            environment=payload.environment.model_dump(),
            weights=payload.weights.model_dump(),
            n_decoys=payload.n_decoys,
            n_random=payload.n_random,
            percent_min=payload.percent_min,
            percent_max=payload.percent_max,
            percent_step=payload.percent_step,
            seed=payload.seed,
        )
        _log_usage_event(request, "optimal-penetration", payload, "ok", time.perf_counter() - started)
        return result
    except Exception as exc:
        _log_usage_event(request, "optimal-penetration", payload, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/pdb")
def pdb(payload: PdbRequest, request: Request) -> Dict[str, Any]:
    started = time.perf_counter()
    try:
        result = {"sequence": payload.sequence.upper(), "pdb": generate_helix_pdb(payload.sequence)}
        _log_usage_event(request, "pdb", payload, "ok", time.perf_counter() - started)
        return result
    except Exception as exc:
        _log_usage_event(request, "pdb", payload, "error", time.perf_counter() - started, {"error": str(exc)})
        raise HTTPException(status_code=400, detail=str(exc)) from exc
