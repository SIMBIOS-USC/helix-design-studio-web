from __future__ import annotations

import math
import os
import statistics
import sys
import tempfile
from pathlib import Path
from typing import Any, Callable, Dict, Generator, List, Optional

import numpy as np

def _resolve_code_dir() -> Path:
    env_override = os.environ.get("QFOLD_CODE_DIR")
    candidate = (Path(env_override).expanduser().resolve() if env_override
                 else Path(__file__).resolve().parents[1] / "Code")
    if not (candidate / "core" / "hamiltonian_builder.py").is_file():
        raise FileNotFoundError(f"Could not locate the Helix Design Studio runtime at {candidate}.")
    return candidate


CODE_DIR = _resolve_code_dir()
if str(CODE_DIR) not in sys.path:
    sys.path.insert(0, str(CODE_DIR))

from core.hamiltonian_builder import HamiltonianBuilder  # noqa: E402

DEFAULT_ALPHABET_16 = ["A", "P", "D", "E", "Q", "L", "M", "F", "K", "S", "T", "W", "Y", "V", "H", "G"]
DEFAULT_WEIGHTS = {
    "lambda_env": 6.0,
    "lambda_charge": 3.0,
    "lambda_local": 1.5,
    "lambda_pairwise": 0.1,
    "lambda_helix_pairs": 1.5,
    "lambda_electrostatic": 4.0,
}
DEFAULT_CACHE_DIR = os.environ.get(
    "HELIX_CACHE_DIR", str(Path(tempfile.gettempdir()) / "helix-design-studio-cache")
)

ENVIRONMENT_PRESETS: Dict[str, Dict[str, Any]] = {
    "homogeneous_polar": {
        "label": "Homogeneous polar medium",
        "wheel_phase_deg": 90.0,
        "wheel_halfwidth_deg": -1.0,
        "membrane_charge": "neu",
    },
    "homogeneous_apolar": {
        "label": "Homogeneous apolar medium",
        "wheel_phase_deg": 90.0,
        "wheel_halfwidth_deg": 180.0,
        "membrane_charge": "neu",
    },
    "interfacial_neg": {
        "label": "Interfacial, anionic membrane",
        "wheel_phase_deg": 90.0,
        "wheel_halfwidth_deg": 90.0,
        "membrane_charge": "neg",
    },
    "interfacial_neu": {
        "label": "Interfacial, neutral membrane",
        "wheel_phase_deg": 90.0,
        "wheel_halfwidth_deg": 90.0,
        "membrane_charge": "neu",
    },
    "interfacial_pos": {
        "label": "Interfacial, cationic membrane",
        "wheel_phase_deg": 90.0,
        "wheel_halfwidth_deg": 90.0,
        "membrane_charge": "pos",
    },
}

HYDROPHOBICITY_FP = {
    "D": -0.77, "E": -0.64, "K": -0.99, "R": -1.01, "H": 0.13,
    "G": 0.00, "A": 0.31, "V": 1.22, "L": 1.70, "I": 1.80,
    "P": 0.72, "M": 1.23, "F": 1.79, "W": 2.25, "Y": 0.96,
    "T": 0.26, "S": -0.04, "C": 1.54, "N": -0.60, "Q": -0.22,
}

RESIDUE_COLORS = {
    "hydrophobic": "#b05a1b",
    "polar": "#2b8a78",
    "positive": "#e03131",
    "negative": "#265dcb",
    "special": "#7a5cfa",
}

HYDROPHOBIC = {"A", "V", "L", "I", "M", "F", "W", "Y", "C", "P"}
POLAR = {"Q", "N", "S", "T", "G"}
POSITIVE = {"K", "R"}
NEGATIVE = {"D", "E"}
SPECIAL = {"H", "X"}

THREE_LETTER_CODES = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}


def bits_per_position(n_aa: int) -> int:
    return max(1, math.ceil(math.log2(n_aa)))


def normalize_residues(residues: Optional[List[str] | str]) -> List[str]:
    if residues is None:
        return list(DEFAULT_ALPHABET_16)
    if isinstance(residues, str):
        if "," in residues:
            toks = [t.strip().upper() for t in residues.split(",") if t.strip()]
        else:
            toks = [c.strip().upper() for c in residues if c.strip()]
    else:
        toks = [str(t).strip().upper() for t in residues if str(t).strip()]
    seen: List[str] = []
    for aa in toks:
        if aa not in seen:
            seen.append(aa)
    return seen or list(DEFAULT_ALPHABET_16)


def normalize_weights(weights: Optional[Dict[str, float]]) -> Dict[str, float]:
    merged = dict(DEFAULT_WEIGHTS)
    if weights:
        merged.update({k: float(v) for k, v in weights.items() if v is not None})
    return merged


def resolve_environment(environment: Optional[Dict[str, Any]]) -> Dict[str, Any]:
    environment = environment or {}
    preset = environment.get("preset", "interfacial_neg")
    base = dict(ENVIRONMENT_PRESETS.get(preset, ENVIRONMENT_PRESETS["interfacial_neg"]))
    for key in ("wheel_phase_deg", "wheel_halfwidth_deg", "membrane_charge"):
        if environment.get(key) is not None:
            base[key] = environment[key]
    base["preset"] = preset
    return base


def align_interfacial_environment_to_hydrophobic_moment(
    sequence: str,
    environment: Dict[str, Any],
) -> Dict[str, Any]:
    descriptor = hydrophobic_moment_descriptor(sequence)
    env = dict(environment)
    if not str(env.get("preset", "")).startswith("interfacial_"):
        return env

    if descriptor["magnitude"] < 1e-10:
        env["phase_aligned_to_hydrophobic_moment"] = False
        env["hydrophobic_moment_phase_deg"] = 0.0
        return env

    native_phase_deg = descriptor["phase_deg"]
    aligned_phase_deg = (-native_phase_deg) % 360.0
    env["wheel_phase_deg"] = aligned_phase_deg
    env["phase_aligned_to_hydrophobic_moment"] = True
    env["hydrophobic_moment_phase_deg"] = native_phase_deg
    return env


def hydrophobic_moment_descriptor(sequence: str) -> Dict[str, float]:
    mx = 0.0
    my = 0.0
    cleaned = sequence.strip().upper()
    for idx, aa in enumerate(cleaned, start=1):
        hydro = HYDROPHOBICITY_FP.get(aa, 0.0)
        theta_deg = (idx - 1) * 100.0
        theta = math.radians(theta_deg)
        mx += hydro * math.cos(theta)
        my += hydro * math.sin(theta)
    magnitude = math.hypot(mx, my)
    phase_deg = math.degrees(math.atan2(my, mx)) if magnitude >= 1e-10 else 0.0
    return {
        "x": float(mx),
        "y": float(my),
        "magnitude": float(magnitude),
        "magnitude_per_residue": float(magnitude / len(cleaned)) if cleaned else 0.0,
        "phase_deg": float(phase_deg),
    }


def penetration_percent_to_halfwidth(percent: float) -> float:
    clamped = max(0.0, min(100.0, float(percent)))
    return 40.0 + (clamped / 100.0) * 110.0


def penetration_percent_to_environment(percent: int, interface_preset: str) -> Dict[str, Any]:
    clamped = max(0, min(100, int(percent)))
    if clamped <= 0:
        env = dict(ENVIRONMENT_PRESETS["homogeneous_polar"])
        env["preset"] = "homogeneous_polar"
        env["penetration_regime"] = "polar_endpoint"
        return env
    if clamped >= 100:
        env = dict(ENVIRONMENT_PRESETS["homogeneous_apolar"])
        env["preset"] = "homogeneous_apolar"
        env["penetration_regime"] = "apolar_endpoint"
        return env
    env = dict(ENVIRONMENT_PRESETS.get(interface_preset, ENVIRONMENT_PRESETS["interfacial_neg"]))
    env["preset"] = interface_preset
    env["wheel_halfwidth_deg"] = penetration_percent_to_halfwidth(clamped)
    env["penetration_regime"] = "interfacial"
    return env


def build_builder(
    length: int,
    amino_acids: List[str],
    environment: Dict[str, Any],
    n_decoys: int = 2000,
    zscore_seed: int = 42,
    stats_cache_dir: str = DEFAULT_CACHE_DIR,
) -> HamiltonianBuilder:
    bpp = bits_per_position(len(amino_acids))
    return HamiltonianBuilder(
        L=length,
        amino_acids=amino_acids,
        bits_per_pos=bpp,
        n_qubits=length * bpp,
        n_decoys=n_decoys,
        zscore_seed=zscore_seed,
        stats_cache_dir=stats_cache_dir,
        membrane_mode="wheel",
        wheel_phase_deg=float(environment["wheel_phase_deg"]),
        wheel_halfwidth_deg=float(environment["wheel_halfwidth_deg"]),
        membrane_charge=str(environment["membrane_charge"]),
        max_interaction_dist=1,
    )


def seq_to_codes(sequence: str, amino_acids: List[str]) -> List[int]:
    aa_to_idx = {aa: i for i, aa in enumerate(amino_acids)}
    try:
        return [aa_to_idx[aa] for aa in sequence]
    except KeyError as exc:
        raise ValueError(f"Residue {exc.args[0]} is not present in the selected alphabet") from exc


def codes_to_seq(codes: List[int], amino_acids: List[str]) -> str:
    return "".join(amino_acids[idx] for idx in codes)


def compute_breakdown(builder: HamiltonianBuilder, codes: List[int], weights: Dict[str, float]) -> Dict[str, float]:
    raw = {
        "local": builder._raw_helix_local(codes),
        "env": builder._raw_env_pol(codes),
        "charge": builder._raw_env_chg(codes),
        "pairwise": builder._raw_pw_int(codes),
        "electrostatic": builder._raw_electrostatic(codes),
        "neigh": builder._raw_helix_neigh(codes),
    }
    z = {
        "local": builder._apply_z_score("helix_local", raw["local"]),
        "env": builder._apply_z_score("env_pol", raw["env"]),
        "charge": builder._apply_z_score("env_chg", raw["charge"]),
        "pairwise": builder._apply_z_score("pw_int", raw["pairwise"]),
        "electrostatic": builder._apply_z_score("electrostatic", raw["electrostatic"]),
        "neigh": builder._apply_z_score("helix_neigh", raw["neigh"]),
    }
    weighted = {
        "local": weights["lambda_local"] * z["local"],
        "env": weights["lambda_env"] * z["env"],
        "charge": weights["lambda_charge"] * z["charge"],
        "pairwise": weights["lambda_pairwise"] * z["pairwise"],
        "electrostatic": weights["lambda_electrostatic"] * z["electrostatic"],
        "neigh": weights["lambda_helix_pairs"] * z["neigh"],
    }
    weighted["total"] = float(sum(weighted.values()))
    return weighted


def residue_family(aa: str) -> str:
    if aa in POSITIVE:
        return "positive"
    if aa in NEGATIVE:
        return "negative"
    if aa in POLAR:
        return "polar"
    if aa in HYDROPHOBIC:
        return "hydrophobic"
    return "special"


def wheel_data(sequence: str, environment: Dict[str, Any]) -> List[Dict[str, Any]]:
    phase = float(environment["wheel_phase_deg"])
    half_width = float(environment["wheel_halfwidth_deg"])
    data = []
    for idx, aa in enumerate(sequence, start=1):
        theta = ((idx - 1) * 100.0 + phase) % 360.0
        signed = theta if theta <= 180.0 else theta - 360.0
        env = "membrane" if abs(signed) <= half_width else "water"
        rad = math.radians(theta - 90.0)
        data.append(
            {
                "index": idx,
                "residue": aa,
                "angle_deg": theta,
                "signed_angle_deg": signed,
                "x": round(math.cos(rad), 6),
                "y": round(math.sin(rad), 6),
                "environment": env,
                "family": residue_family(aa),
                "color": RESIDUE_COLORS[residue_family(aa)],
            }
        )
    return data


def model_net_charge(sequence: str) -> int:
    charge_map = {"R": 1, "K": 1, "D": -1, "E": -1}
    return int(sum(charge_map.get(aa, 0) for aa in sequence))


def neutral_ph_net_charge(sequence: str) -> int:
    charge_map = {"R": 1, "K": 1, "D": -1, "E": -1}
    return int(sum(charge_map.get(aa, 0) for aa in sequence))


def random_reference(
    builder: HamiltonianBuilder,
    amino_acids: List[str],
    weights: Dict[str, float],
    n_random: int,
    seed: int,
) -> Dict[str, float]:
    rng = np.random.default_rng(seed)
    energies = []
    for _ in range(n_random):
        codes = rng.integers(0, len(amino_acids), size=builder.L).tolist()
        energies.append(compute_breakdown(builder, codes, weights)["total"])
    mean = statistics.fmean(energies)
    std = statistics.pstdev(energies) or 1.0
    return {"mean": float(mean), "std": float(std), "energies": energies, "seed": seed}


def score_sequence_with_builder(
    sequence: str,
    builder: HamiltonianBuilder,
    environment: Dict[str, Any],
    amino_acids: List[str],
    weights: Dict[str, float],
    ref: Dict[str, float],
    n_random: int,
) -> Dict[str, Any]:
    sequence = sequence.strip().upper()
    codes = seq_to_codes(sequence, amino_acids)
    breakdown = compute_breakdown(builder, codes, weights)
    mu_descriptor = hydrophobic_moment_descriptor(sequence)
    energy = breakdown["total"]
    z_score = (energy - ref["mean"]) / ref["std"]
    percentile = 100.0 * sum(e <= energy for e in ref["energies"]) / len(ref["energies"])
    return {
        "sequence": sequence,
        "length": len(sequence),
        "environment": environment,
        "weights": weights,
        "energy": energy,
        "z_score": float(z_score),
        "percentile_lower_is_better": float(percentile),
        "random_reference": {
            "mean": ref["mean"], "std": ref["std"], "n": n_random, "seed": ref["seed"],
        },
        "calibration": {
            "n_decoys": builder.kwargs.get("n_decoys", 2000),
            "seed": builder.kwargs.get("zscore_seed", 42),
            "cache_version": builder.stats_cache_version,
        },
        "orientation_mode": (
            "hydrophobic_moment" if "phase_aligned_to_hydrophobic_moment" in environment else "fixed"
        ),
        "breakdown": breakdown,
        "observables": {
            "abs_mu_h": mu_descriptor["magnitude"],
            "abs_mu_h_per_residue": mu_descriptor["magnitude_per_residue"],
            "mu_h_phase_deg": mu_descriptor["phase_deg"],
        },
        "wheel": wheel_data(sequence, environment),
        "model_net_charge": model_net_charge(sequence),
        "neutral_ph_net_charge": neutral_ph_net_charge(sequence),
        "alphabet": amino_acids,
    }


def hamming_distance_codes(a: List[int], b: List[int]) -> int:
    return sum(x != y for x, y in zip(a, b))


def diversity_penalty(
    codes: List[int],
    accepted_codes: List[List[int]],
    min_distance: int,
    penalty_weight: float,
) -> float:
    if not accepted_codes:
        return 0.0
    penalty = 0.0
    for prior in accepted_codes:
        distance = hamming_distance_codes(codes, prior)
        deficit = max(0, min_distance - distance)
        if deficit:
            penalty += penalty_weight * deficit
    return penalty


def score_sequence(
    sequence: str,
    environment: Optional[Dict[str, Any]] = None,
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    n_random: int = 500,
    seed: int = 42,
) -> Dict[str, Any]:
    sequence = sequence.strip().upper()
    amino_acids = normalize_residues(residues)
    env = resolve_environment(environment)
    env = align_interfacial_environment_to_hydrophobic_moment(sequence, env)
    w = normalize_weights(weights)
    builder = build_builder(len(sequence), amino_acids, env, n_decoys=n_decoys, zscore_seed=seed)
    ref = random_reference(builder, amino_acids, w, n_random=n_random, seed=seed + 17)
    return score_sequence_with_builder(sequence, builder, env, amino_acids, w, ref, n_random)


def score_sequence_fixed_environment(
    sequence: str,
    environment: Optional[Dict[str, Any]] = None,
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    n_random: int = 500,
    seed: int = 42,
) -> Dict[str, Any]:
    sequence = sequence.strip().upper()
    amino_acids = normalize_residues(residues)
    env = resolve_environment(environment)
    w = normalize_weights(weights)
    builder = build_builder(len(sequence), amino_acids, env, n_decoys=n_decoys, zscore_seed=seed)
    ref = random_reference(builder, amino_acids, w, n_random=n_random, seed=seed + 17)
    return score_sequence_with_builder(sequence, builder, env, amino_acids, w, ref, n_random)


def simulated_annealing(
    builder: HamiltonianBuilder,
    amino_acids: List[str],
    weights: Dict[str, float],
    steps: int,
    restarts: int,
    seed: int,
    objective_fn: Optional[Callable[[List[int]], float]] = None,
) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    n_aa = len(amino_acids)
    best_codes = None
    best_objective = float("inf")
    best_energy = float("inf")

    for _ in range(restarts):
        current = rng.integers(0, n_aa, size=builder.L).tolist()
        current_energy = compute_breakdown(builder, current, weights)["total"]
        current_objective = objective_fn(current) if objective_fn else current_energy
        local_best = list(current)
        local_best_energy = current_energy
        local_best_objective = current_objective

        for step in range(steps):
            temperature = max(1e-4, 2.0 * (1.0 - step / steps))
            proposal = list(current)
            pos = int(rng.integers(0, builder.L))
            old_code = proposal[pos]
            new_code = int(rng.integers(0, n_aa - 1))
            if new_code >= old_code:
                new_code += 1
            proposal[pos] = new_code

            proposal_energy = compute_breakdown(builder, proposal, weights)["total"]
            proposal_objective = objective_fn(proposal) if objective_fn else proposal_energy
            delta = proposal_objective - current_objective
            if delta <= 0 or rng.random() < math.exp(-delta / temperature):
                current = proposal
                current_energy = proposal_energy
                current_objective = proposal_objective
                if current_objective < local_best_objective:
                    local_best = list(current)
                    local_best_energy = current_energy
                    local_best_objective = current_objective

        if local_best_objective < best_objective:
            best_codes = list(local_best)
            best_energy = local_best_energy
            best_objective = local_best_objective

    assert best_codes is not None
    sequence = codes_to_seq(best_codes, amino_acids)
    return {
        "sequence": sequence,
        "energy": best_energy,
        "objective": best_objective,
        "breakdown": compute_breakdown(builder, best_codes, weights),
        "codes": list(best_codes),
    }


def _iter_specificity_search(
    target_builder: HamiltonianBuilder,
    off_builders: List[HamiltonianBuilder],
    amino_acids: List[str],
    weights: Dict[str, float],
    lambda_balance: float,
    steps: int,
    restarts: int,
    seed: int,
    accepted_codes: Optional[List[List[int]]] = None,
    min_distance: int = 0,
    penalty_weight: float = 0.0,
    progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
) -> Dict[str, Any]:
    rng = np.random.default_rng(seed)
    n_aa = len(amino_acids)
    accepted_codes = accepted_codes or []
    best_codes = None
    best_objective = float("inf")
    best_penalized_objective = float("inf")
    best_target_energy = float("inf")
    best_off_energies: List[float] = []
    checkpoint = max(250, min(2000, steps // 30))

    for restart_idx in range(restarts):
        current = rng.integers(0, n_aa, size=target_builder.L).tolist()
        current_target = compute_breakdown(target_builder, current, weights)["total"]
        current_off = [compute_breakdown(builder, current, weights)["total"] for builder in off_builders]
        current_raw_objective = current_target + (lambda_balance - 1.0) * min(current_off)
        current_objective = current_raw_objective + diversity_penalty(
            current,
            accepted_codes,
            min_distance=min_distance,
            penalty_weight=penalty_weight,
        )
        local_best = list(current)
        local_best_target = current_target
        local_best_off = list(current_off)
        local_best_objective = current_raw_objective
        local_best_penalized_objective = current_objective

        progress_payload = {
            "phase": "annealing",
            "restart": restart_idx + 1,
            "restarts": restarts,
            "step": 0,
            "steps": steps,
            "progress_fraction": restart_idx / max(1, restarts),
            "latest_sequence": codes_to_seq(current, amino_acids),
            "best_sequence": codes_to_seq(local_best, amino_acids),
            "best_objective": local_best_objective,
        }
        if progress_callback:
            progress_callback(progress_payload)
        yield {"type": "progress", "payload": progress_payload}

        for step_idx in range(steps):
            temperature = max(1e-4, 2.0 * (1.0 - step_idx / steps))
            proposal = list(current)
            pos = int(rng.integers(0, target_builder.L))
            old_code = proposal[pos]
            new_code = int(rng.integers(0, n_aa - 1))
            if new_code >= old_code:
                new_code += 1
            proposal[pos] = new_code

            proposal_target = compute_breakdown(target_builder, proposal, weights)["total"]
            proposal_off = [compute_breakdown(builder, proposal, weights)["total"] for builder in off_builders]
            proposal_raw_objective = proposal_target + (lambda_balance - 1.0) * min(proposal_off)
            proposal_objective = proposal_raw_objective + diversity_penalty(
                proposal,
                accepted_codes,
                min_distance=min_distance,
                penalty_weight=penalty_weight,
            )
            delta = proposal_objective - current_objective
            if delta <= 0 or rng.random() < math.exp(-delta / temperature):
                current = proposal
                current_target = proposal_target
                current_off = proposal_off
                current_objective = proposal_objective
                current_raw_objective = proposal_raw_objective
                if current_objective < local_best_penalized_objective:
                    local_best = list(current)
                    local_best_target = current_target
                    local_best_off = list(current_off)
                    local_best_objective = current_raw_objective
                    local_best_penalized_objective = current_objective

            if (step_idx + 1) % checkpoint == 0 or step_idx + 1 == steps:
                overall = (restart_idx + (step_idx + 1) / max(1, steps)) / max(1, restarts)
                best_sequence_codes = (
                    best_codes
                    if best_codes is not None and best_penalized_objective <= local_best_penalized_objective
                    else local_best
                )
                best_sequence = codes_to_seq(best_sequence_codes, amino_acids)
                best_obj = min(best_objective, local_best_objective)
                progress_payload = {
                    "phase": "annealing",
                    "restart": restart_idx + 1,
                    "restarts": restarts,
                    "step": step_idx + 1,
                    "steps": steps,
                    "progress_fraction": overall,
                    "latest_sequence": codes_to_seq(current, amino_acids),
                    "best_sequence": best_sequence,
                    "best_objective": best_obj,
                }
                if progress_callback:
                    progress_callback(progress_payload)
                yield {"type": "progress", "payload": progress_payload}

        if local_best_penalized_objective < best_penalized_objective:
            best_codes = list(local_best)
            best_objective = local_best_objective
            best_penalized_objective = local_best_penalized_objective
            best_target_energy = local_best_target
            best_off_energies = list(local_best_off)

    assert best_codes is not None
    yield {
        "type": "result",
        "payload": {
            "sequence": codes_to_seq(best_codes, amino_acids),
            "codes": list(best_codes),
            "objective": best_objective,
            "penalized_objective": best_penalized_objective,
            "best_target_energy": best_target_energy,
            "best_off_energies": best_off_energies,
        },
    }


def _run_specificity_search(
    target_builder: HamiltonianBuilder,
    off_builders: List[HamiltonianBuilder],
    amino_acids: List[str],
    weights: Dict[str, float],
    lambda_balance: float,
    steps: int,
    restarts: int,
    seed: int,
    accepted_codes: Optional[List[List[int]]] = None,
    min_distance: int = 0,
    penalty_weight: float = 0.0,
    progress_callback: Optional[Callable[[Dict[str, Any]], None]] = None,
) -> Dict[str, Any]:
    final_payload: Optional[Dict[str, Any]] = None
    for update in _iter_specificity_search(
        target_builder=target_builder,
        off_builders=off_builders,
        amino_acids=amino_acids,
        weights=weights,
        lambda_balance=lambda_balance,
        steps=steps,
        restarts=restarts,
        seed=seed,
        accepted_codes=accepted_codes,
        min_distance=min_distance,
        penalty_weight=penalty_weight,
        progress_callback=progress_callback,
    ):
        if update["type"] == "result":
            final_payload = update["payload"]
    if final_payload is None:
        raise RuntimeError("Specificity search did not produce a result")
    return final_payload


def _build_specificity_candidate_payload(
    sequence: str,
    *,
    length: int,
    amino_acids: List[str],
    target_env: Dict[str, Any],
    off_target_envs: List[Dict[str, Any]],
    weights: Dict[str, float],
    lambda_balance: float,
    steps: int,
    restarts: int,
    seed: int,
    objective: float,
    best_target_energy: float,
    best_off_energies: List[float],
    target_builder: HamiltonianBuilder,
    off_builders: List[HamiltonianBuilder],
) -> Dict[str, Any]:
    # Use the search Hamiltonians so reported scores retain the same
    # calibration and orientation.
    n_random = 3000
    calibration_seed = int(target_builder.kwargs.get("zscore_seed", 42))
    target_ref = random_reference(
        target_builder, amino_acids, weights, n_random, calibration_seed + 17
    )
    target_score = score_sequence_with_builder(
        sequence, target_builder, target_env, amino_acids, weights, target_ref, n_random
    )
    off_target_scores = []
    for idx, (env, builder) in enumerate(zip(off_target_envs, off_builders)):
        ref = random_reference(builder, amino_acids, weights, n_random, calibration_seed + 18 + idx)
        off_target_scores.append(
            score_sequence_with_builder(
                sequence, builder, env, amino_acids, weights, ref, n_random
            )
        )

    worst_off_target = min(off_target_scores, key=lambda item: item["energy"])
    specificity_design = {
        "objective": objective,
        "objective_formula": "H_target + (lambda - 1) * min(H_offtargets)",
        "objective_formula_expanded": "lambda * H_target + (1 - lambda) * (H_target - min(H_offtargets))",
        "lambda_balance": lambda_balance,
        "target_preset": target_env["preset"],
        "off_target_presets": [env["preset"] for env in off_target_envs],
        "worst_off_target_preset": worst_off_target["environment"]["preset"],
        "energy_margin": worst_off_target["energy"] - target_score["energy"],
        "z_score_margin": worst_off_target["z_score"] - target_score["z_score"],
        "steps": steps,
        "restarts": restarts,
        "seed": seed,
        "best_target_energy": best_target_energy,
        "best_off_target_energy": min(best_off_energies) if best_off_energies else None,
    }
    return {
        "sequence": sequence,
        "length": length,
        "alphabet": amino_acids,
        "target_environment": target_score,
        "off_target_environments": off_target_scores,
        "specificity_design": specificity_design,
    }


def design_sequence(
    length: int,
    environment: Optional[Dict[str, Any]] = None,
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    steps: int = 12000,
    restarts: int = 8,
    seed: int = 42,
) -> Dict[str, Any]:
    amino_acids = normalize_residues(residues)
    env = resolve_environment(environment)
    w = normalize_weights(weights)
    builder = build_builder(length, amino_acids, env, n_decoys=n_decoys, zscore_seed=seed)
    ref = random_reference(builder, amino_acids, w, n_random=500, seed=seed + 17)
    result = simulated_annealing(builder, amino_acids, w, steps=steps, restarts=restarts, seed=seed)
    score = score_sequence_with_builder(result["sequence"], builder, env, amino_acids, w, ref, 500)
    score["design_settings"] = {"steps": steps, "restarts": restarts, "seed": seed}
    return score


def design_family(
    length: int,
    environment: Optional[Dict[str, Any]] = None,
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    steps: int = 10000,
    restarts: int = 6,
    family_size: int = 8,
    oversample_factor: int = 3,
    seed: int = 42,
) -> Dict[str, Any]:
    final_update = None
    for update in design_family_stream(
        length=length,
        environment=environment,
        residues=residues,
        weights=weights,
        n_decoys=n_decoys,
        steps=steps,
        restarts=restarts,
        family_size=family_size,
        oversample_factor=oversample_factor,
        seed=seed,
    ):
        if update["type"] == "final":
            final_update = update["payload"]
    if final_update is None:
        raise RuntimeError("Family generation did not produce a final result")
    return final_update


def design_family_stream(
    length: int,
    environment: Optional[Dict[str, Any]] = None,
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    steps: int = 10000,
    restarts: int = 6,
    family_size: int = 8,
    oversample_factor: int = 3,
    seed: int = 42,
) -> Generator[Dict[str, Any], None, None]:
    amino_acids = normalize_residues(residues)
    env = resolve_environment(environment)
    w = normalize_weights(weights)
    builder = build_builder(length, amino_acids, env, n_decoys=n_decoys, zscore_seed=seed)
    n_random = 500
    ref = random_reference(builder, amino_acids, w, n_random=n_random, seed=seed + 17)

    unique_results: Dict[str, Dict[str, Any]] = {}
    accepted_codes: List[List[int]] = []
    total_runs = max(family_size * max(oversample_factor, 4), family_size + 6)
    min_distance = max(2, math.ceil(length * 0.15))
    penalty_weight = 6.0

    yield {
        "type": "meta",
        "payload": {
            "target": family_size,
            "total_runs_budget": total_runs,
            "min_distance": min_distance,
            "environment": env,
        },
    }

    for run_idx in range(total_runs):
        objective_fn = None
        if accepted_codes:
            objective_fn = lambda codes, prior=accepted_codes: (
                compute_breakdown(builder, codes, w)["total"]
                + diversity_penalty(codes, prior, min_distance=min_distance, penalty_weight=penalty_weight)
            )
        candidate = simulated_annealing(
            builder,
            amino_acids,
            w,
            steps=steps,
            restarts=restarts,
            seed=seed + run_idx,
            objective_fn=objective_fn,
        )
        sequence = candidate["sequence"]
        yield {
            "type": "progress",
            "payload": {
                "attempt": run_idx + 1,
                "accepted": len(unique_results),
                "target": family_size,
                "latest_sequence": sequence,
            },
        }
        if sequence in unique_results:
            continue
        scored = score_sequence_with_builder(sequence, builder, env, amino_acids, w, ref, n_random)
        min_hamming = None
        if accepted_codes:
            min_hamming = min(hamming_distance_codes(candidate["codes"], prior) for prior in accepted_codes)
        unique_results[sequence] = scored
        accepted_codes.append(candidate["codes"])
        yield {
            "type": "candidate",
            "payload": {
                "attempt": run_idx + 1,
                "accepted": len(unique_results),
                "target": family_size,
                "min_hamming_to_previous": min_hamming,
                "member": scored,
            },
        }
        if len(unique_results) >= family_size:
            break

    ranked = sorted(unique_results.values(), key=lambda item: item["energy"])[:family_size]
    payload = {
        "environment": env,
        "weights": w,
        "alphabet": amino_acids,
        "family_size_requested": family_size,
        "family_size_returned": len(ranked),
        "design_settings": {
            "steps": steps,
            "restarts": restarts,
            "seed": seed,
            "oversample_factor": oversample_factor,
            "min_distance": min_distance,
            "penalty_weight": penalty_weight,
        },
        "members": ranked,
        "status_message": (
            f"Generated {len(ranked)} unique candidates"
            if len(ranked) >= family_size
            else f"Generated {len(ranked)} unique candidates out of the requested {family_size}"
        ),
    }
    yield {"type": "final", "payload": payload}


def compare_environments(
    sequence: str,
    environment_a: Dict[str, Any],
    environment_b: Dict[str, Any],
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    seed: int = 42,
) -> Dict[str, Any]:
    amino_acids = normalize_residues(residues)
    w = normalize_weights(weights)
    env_a = resolve_environment(environment_a)
    env_b = resolve_environment(environment_b)
    scored_a = score_sequence(sequence, env_a, amino_acids, w, n_decoys=n_decoys, seed=seed)
    scored_b = score_sequence(sequence, env_b, amino_acids, w, n_decoys=n_decoys, seed=seed + 1)
    return {
        "sequence": sequence.upper(),
        "environment_a": scored_a,
        "environment_b": scored_b,
        "energy_gap": abs(scored_a["energy"] - scored_b["energy"]),
        "z_score_gap": abs(scored_a["z_score"] - scored_b["z_score"]),
        "joint_energy": scored_a["energy"] + scored_b["energy"],
    }


def compare_environments_fixed(
    sequence: str,
    environment_a: Dict[str, Any],
    environment_b: Dict[str, Any],
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    n_random: int = 500,
    seed: int = 42,
) -> Dict[str, Any]:
    amino_acids = normalize_residues(residues)
    w = normalize_weights(weights)
    env_a = resolve_environment(environment_a)
    env_b = resolve_environment(environment_b)
    builder_a = build_builder(len(sequence), amino_acids, env_a, n_decoys=n_decoys, zscore_seed=seed)
    builder_b = build_builder(len(sequence), amino_acids, env_b, n_decoys=n_decoys, zscore_seed=seed + 1)
    ref_a = random_reference(builder_a, amino_acids, w, n_random=n_random, seed=seed + 17)
    ref_b = random_reference(builder_b, amino_acids, w, n_random=n_random, seed=seed + 18)
    scored_a = score_sequence_with_builder(sequence, builder_a, env_a, amino_acids, w, ref_a, n_random)
    scored_b = score_sequence_with_builder(sequence, builder_b, env_b, amino_acids, w, ref_b, n_random)
    return {
        "sequence": sequence.upper(),
        "environment_a": scored_a,
        "environment_b": scored_b,
        "energy_gap": abs(scored_a["energy"] - scored_b["energy"]),
        "z_score_gap": abs(scored_a["z_score"] - scored_b["z_score"]),
        "joint_energy": scored_a["energy"] + scored_b["energy"],
    }


def optimize_penetration(
    sequence: str,
    environment: Dict[str, Any],
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    n_random: int = 500,
    percent_min: int = 0,
    percent_max: int = 100,
    percent_step: int = 5,
    seed: int = 42,
) -> Dict[str, Any]:
    amino_acids = normalize_residues(residues)
    w = normalize_weights(weights)
    env_base = resolve_environment(environment)
    if not str(env_base["preset"]).startswith("interfacial_"):
        raise ValueError("Optimal penetration is only defined for interfacial presets.")
    if percent_step <= 0:
        raise ValueError("percent_step must be positive.")
    if percent_min > percent_max:
        raise ValueError("percent_min must be smaller than or equal to percent_max.")

    profile: List[Dict[str, Any]] = []
    best_result: Optional[Dict[str, Any]] = None
    best_percent: Optional[int] = None

    for idx, percent in enumerate(range(percent_min, percent_max + 1, percent_step)):
        env = penetration_percent_to_environment(percent, str(env_base["preset"]))
        scored = score_sequence(
            sequence=sequence,
            environment=env,
            residues=amino_acids,
            weights=w,
            n_decoys=n_decoys,
            n_random=n_random,
            seed=seed + idx,
        )
        profile.append(
            {
                "penetration_percent": percent,
                "wheel_halfwidth_deg": env.get("wheel_halfwidth_deg"),
                "environment_preset": env["preset"],
                "penetration_regime": env.get("penetration_regime", "interfacial"),
                "energy": scored["energy"],
                "z_score": scored["z_score"],
                "percentile_lower_is_better": scored["percentile_lower_is_better"],
            }
        )
        if best_result is None or scored["z_score"] < best_result["z_score"]:
            best_result = scored
            best_percent = percent

    assert best_result is not None and best_percent is not None
    return {
        "sequence": sequence.strip().upper(),
        "interface_preset": env_base["preset"],
        "optimization_target": "minimum_z_score",
        "best_penetration_percent": best_percent,
        "best_score": best_result,
        "profile": profile,
        "scan_settings": {
            "percent_min": percent_min,
            "percent_max": percent_max,
            "percent_step": percent_step,
            "n_random": n_random,
            "n_decoys": n_decoys,
            "seed": seed,
        },
    }


def cross_environment_design(
    length: int,
    environment_a: Dict[str, Any],
    environment_b: Dict[str, Any],
    lambda_gap: float,
    residues: Optional[List[str] | str] = None,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 2000,
    steps: int = 16000,
    restarts: int = 10,
    seed: int = 42,
) -> Dict[str, Any]:
    amino_acids = normalize_residues(residues)
    w = normalize_weights(weights)
    env_a = resolve_environment(environment_a)
    env_b = resolve_environment(environment_b)
    builder_a = build_builder(length, amino_acids, env_a, n_decoys=n_decoys, zscore_seed=seed)
    builder_b = build_builder(length, amino_acids, env_b, n_decoys=n_decoys, zscore_seed=seed)
    rng = np.random.default_rng(seed)
    n_aa = len(amino_acids)
    best_codes = None
    best_obj = float("inf")
    best_a = 0.0
    best_b = 0.0

    for _ in range(restarts):
        current = rng.integers(0, n_aa, size=length).tolist()
        current_a = compute_breakdown(builder_a, current, w)["total"]
        current_b = compute_breakdown(builder_b, current, w)["total"]
        current_obj = current_a + current_b + lambda_gap * abs(current_a - current_b)
        local_best = list(current)
        local_obj = current_obj
        local_a = current_a
        local_b = current_b

        for step in range(steps):
            temperature = max(1e-4, 2.0 * (1.0 - step / steps))
            proposal = list(current)
            pos = int(rng.integers(0, length))
            old_code = proposal[pos]
            new_code = int(rng.integers(0, n_aa - 1))
            if new_code >= old_code:
                new_code += 1
            proposal[pos] = new_code
            prop_a = compute_breakdown(builder_a, proposal, w)["total"]
            prop_b = compute_breakdown(builder_b, proposal, w)["total"]
            prop_obj = prop_a + prop_b + lambda_gap * abs(prop_a - prop_b)
            delta = prop_obj - current_obj
            if delta <= 0 or rng.random() < math.exp(-delta / temperature):
                current = proposal
                current_obj = prop_obj
                current_a = prop_a
                current_b = prop_b
                if current_obj < local_obj:
                    local_best = list(current)
                    local_obj = current_obj
                    local_a = current_a
                    local_b = current_b

        if local_obj < best_obj:
            best_codes = list(local_best)
            best_obj = local_obj
            best_a = local_a
            best_b = local_b

    assert best_codes is not None
    sequence = codes_to_seq(best_codes, amino_acids)
    # Reuse both search builders, including their calibration seeds, for display.
    ref_a = random_reference(builder_a, amino_acids, w, n_random=500, seed=seed + 17)
    ref_b = random_reference(builder_b, amino_acids, w, n_random=500, seed=seed + 18)
    scored_a = score_sequence_with_builder(sequence, builder_a, env_a, amino_acids, w, ref_a, 500)
    scored_b = score_sequence_with_builder(sequence, builder_b, env_b, amino_acids, w, ref_b, 500)
    comparison = {
        "sequence": sequence,
        "environment_a": scored_a,
        "environment_b": scored_b,
        "energy_gap": abs(scored_a["energy"] - scored_b["energy"]),
        "z_score_gap": abs(scored_a["z_score"] - scored_b["z_score"]),
        "joint_energy": scored_a["energy"] + scored_b["energy"],
    }
    comparison["cross_design"] = {
        "lambda_gap": lambda_gap,
        "objective": best_obj,
        "energy_a": best_a,
        "energy_b": best_b,
        "steps": steps,
        "restarts": restarts,
        "seed": seed,
    }
    return comparison


def design_specificity(
    length: int,
    target_environment: Dict[str, Any],
    off_target_environments: List[Dict[str, Any]],
    residues: Optional[List[str] | str] = None,
    lambda_balance: float = 0.0,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 10000,
    steps: int = 60000,
    restarts: int = 20,
    num_sequences: int = 1,
    seed: int = 42,
) -> Dict[str, Any]:
    final_update = None
    for update in design_specificity_stream(
        length=length,
        target_environment=target_environment,
        off_target_environments=off_target_environments,
        residues=residues,
        lambda_balance=lambda_balance,
        weights=weights,
        n_decoys=n_decoys,
        steps=steps,
        restarts=restarts,
        num_sequences=num_sequences,
        seed=seed,
    ):
        if update["type"] == "final":
            final_update = update["payload"]
    if final_update is None:
        raise RuntimeError("Specificity design did not produce a final result")
    return final_update


def design_specificity_stream(
    length: int,
    target_environment: Dict[str, Any],
    off_target_environments: List[Dict[str, Any]],
    residues: Optional[List[str] | str] = None,
    lambda_balance: float = 0.0,
    weights: Optional[Dict[str, float]] = None,
    n_decoys: int = 10000,
    steps: int = 60000,
    restarts: int = 20,
    num_sequences: int = 1,
    seed: int = 42,
) -> Generator[Dict[str, Any], None, None]:
    amino_acids = normalize_residues(residues)
    w = normalize_weights(weights)
    target_env = resolve_environment(target_environment)
    off_target_envs = [resolve_environment(env) for env in off_target_environments]
    unique_off_targets: List[Dict[str, Any]] = []
    seen: set[tuple[Any, ...]] = set()
    for env in off_target_envs:
        key = (
            env.get("preset"),
            env.get("wheel_phase_deg"),
            env.get("wheel_halfwidth_deg"),
            env.get("membrane_charge"),
        )
        if key in seen:
            continue
        seen.add(key)
        unique_off_targets.append(env)
    off_target_envs = [
        env for env in unique_off_targets
        if (
            env.get("preset") != target_env.get("preset")
            or float(env.get("wheel_halfwidth_deg", 0.0)) != float(target_env.get("wheel_halfwidth_deg", 0.0))
            or str(env.get("membrane_charge")) != str(target_env.get("membrane_charge"))
        )
    ]
    if not off_target_envs:
        raise ValueError("Please provide at least one off-target environment different from the target.")

    lambda_balance = float(lambda_balance)
    if not (0.0 <= lambda_balance <= 1.0):
        raise ValueError("lambda_balance must lie between 0 and 1.")
    num_sequences = max(1, int(num_sequences))
    total_restart_budget = max(1, int(restarts))
    restarts_per_attempt = max(1, total_restart_budget // num_sequences)
    total_attempts = max(1, total_restart_budget // restarts_per_attempt)

    yield {
        "type": "meta",
        "payload": {
            "target_environment": target_env,
            "off_target_environments": off_target_envs,
            "alphabet_size": len(amino_acids),
            "steps": steps,
            "restarts": restarts_per_attempt,
            "total_restart_budget": total_restart_budget,
            "restarts_per_attempt": restarts_per_attempt,
            "num_sequences": num_sequences,
            "lambda_balance": lambda_balance,
            "phase": "initializing",
            "progress_fraction": 0.0,
        },
    }

    builder_specs = [("target", target_env)] + [(f"off_target_{idx + 1}", env) for idx, env in enumerate(off_target_envs)]
    built: Dict[str, HamiltonianBuilder] = {}
    total_builders = len(builder_specs)
    for idx, (label, env) in enumerate(builder_specs, start=1):
        built[label] = build_builder(length, amino_acids, env, n_decoys=n_decoys, zscore_seed=seed)
        yield {
            "type": "phase",
            "payload": {
                "phase": "calibrating",
                "label": label,
                "completed": idx,
                "total": total_builders,
                "progress_fraction": 0.1 * (idx / max(1, total_builders)),
            },
        }

    target_builder = built["target"]
    off_builders = [built[f"off_target_{idx + 1}"] for idx in range(len(off_target_envs))]
    unique_results: Dict[str, Dict[str, Any]] = {}
    accepted_codes: List[List[int]] = []
    min_distance = max(2, math.ceil(length * 0.15))
    penalty_weight = 6.0

    for attempt_idx in range(total_attempts):
        yield {
            "type": "phase",
            "payload": {
                "phase": "collecting",
                "completed": len(unique_results),
                "total": num_sequences,
                "attempt": attempt_idx + 1,
                "attempts_total": total_attempts,
                "progress_fraction": 0.1 + 0.55 * (attempt_idx / max(1, total_attempts)),
            },
        }

        result: Optional[Dict[str, Any]] = None
        for update in _iter_specificity_search(
            target_builder=target_builder,
            off_builders=off_builders,
            amino_acids=amino_acids,
            weights=w,
            lambda_balance=lambda_balance,
            steps=steps,
            restarts=restarts_per_attempt,
            seed=seed + attempt_idx,
            accepted_codes=accepted_codes,
            min_distance=min_distance,
            penalty_weight=penalty_weight,
        ):
            if update["type"] == "progress":
                inner = update["payload"]
                inner_fraction = float(inner.get("progress_fraction", 0.0))
                overall = 0.1 + 0.55 * ((attempt_idx + inner_fraction) / max(1, total_attempts))
                yield {
                    "type": "progress",
                    "payload": {
                        **inner,
                        "phase": "collecting",
                        "accepted": len(unique_results),
                        "target_count": num_sequences,
                        "attempt": attempt_idx + 1,
                        "attempts_total": total_attempts,
                        "progress_fraction": overall,
                    },
                }
            elif update["type"] == "result":
                result = update["payload"]

        if result is None:
            raise RuntimeError("Specificity search attempt did not return a result")

        sequence = result["sequence"]
        if sequence in unique_results:
            continue

        candidate_payload = _build_specificity_candidate_payload(
            sequence,
            length=length,
            amino_acids=amino_acids,
            target_env=target_env,
            off_target_envs=off_target_envs,
            weights=w,
            lambda_balance=lambda_balance,
            steps=steps,
            restarts=restarts_per_attempt,
            seed=seed + attempt_idx,
            objective=result["objective"],
            best_target_energy=result["best_target_energy"],
            best_off_energies=result["best_off_energies"],
            target_builder=target_builder,
            off_builders=off_builders,
        )
        unique_results[sequence] = candidate_payload
        accepted_codes.append(result["codes"])
        yield {
            "type": "candidate",
            "payload": {
                "accepted": len(unique_results),
                "target": num_sequences,
                "attempt": attempt_idx + 1,
                "sequence": sequence,
                "objective": candidate_payload["specificity_design"]["objective"],
                "z_score_margin": candidate_payload["specificity_design"]["z_score_margin"],
                "member": candidate_payload,
            },
        }
        if len(unique_results) >= num_sequences:
            break

    ranked_candidates = sorted(
        unique_results.values(),
        key=lambda item: (
            round(float(item["specificity_design"]["objective"]), 8),
            -float(item["specificity_design"]["z_score_margin"]),
            float(item["target_environment"]["z_score"]),
            item["sequence"],
        ),
    )
    if not ranked_candidates:
        raise RuntimeError("Specificity design did not recover any candidate sequence.")

    for idx, candidate in enumerate(ranked_candidates, start=1):
        candidate["specificity_design"]["rank"] = idx

    top_candidate = ranked_candidates[0]
    payload = {
        **top_candidate,
        "candidates": ranked_candidates,
        "specificity_collection": {
            "requested": num_sequences,
            "returned": len(ranked_candidates),
            "attempts_budget": total_attempts,
            "total_restart_budget": total_restart_budget,
            "restarts_per_attempt": restarts_per_attempt,
            "lambda_balance": lambda_balance,
            "target_preset": target_env["preset"],
            "off_target_presets": [env["preset"] for env in off_target_envs],
        },
    }
    yield {"type": "final", "payload": payload}
def generate_helix_pdb(sequence: str) -> str:
    sequence = sequence.strip().upper()
    if not sequence:
        raise RuntimeError("Sequence cannot be empty")
    for aa in sequence:
        if aa not in THREE_LETTER_CODES:
            raise RuntimeError(f"Unsupported residue for PDB export: {aa}")

    try:
        import PeptideBuilder
        from Bio.PDB import PDBIO
    except Exception as exc:
        raise RuntimeError(
            "PDB export requires PeptideBuilder and Biopython in the application environment."
        ) from exc

    phi = [-57.0] * max(0, len(sequence) - 1)
    psi = [-47.0] * max(0, len(sequence) - 1)
    structure = PeptideBuilder.make_structure(sequence, phi, psi)

    with tempfile.TemporaryDirectory(prefix="helix_pdb_") as tmpdir:
        out_path = Path(tmpdir) / "helix.pdb"
        io = PDBIO()
        io.set_structure(structure)
        io.save(str(out_path))
        pdb_text = out_path.read_text()
    if pdb_text.strip():
        return pdb_text if pdb_text.endswith("\n") else pdb_text + "\n"
    raise RuntimeError("Failed to generate PDB for the requested sequence.")
