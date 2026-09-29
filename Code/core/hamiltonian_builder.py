"""Classical scoring and calibration used by the Helix Design Studio web app.

The unused Pauli-operator construction backend has been removed from this
web-only distribution. Parameter matrices are indexed by residue identity;
see PARAMETER_PROVENANCE.md for the provisional M3/M4 header reconstruction.
"""

import hashlib
import json
import os
from typing import Any, Dict, List, Tuple

import numpy as np

from data_loaders.energy_matrix_loader import (
    _load_first_neighbors_matrix_file,
    _load_third_neighbors_matrix_file,
    _load_fourth_neighbors_matrix_file,
    _load_energy_matrix_file
)


class HamiltonianBuilder:
    """
    Builds the protein Hamiltonian using intensive property scaling
    and Z-score normalization (Goldstein-Wolynes method).

    Fiel a:
      - Escalado intensivo (Imagen 3 / Sec. Thermodynamic Consistency)
      - Z-score real con N decoys aleatorios por término
      - Momento hidrofóbico con componentes x e y  (|μ_H|² = μx² + μy²)
      - Helix-neighbors sin signo negativo arbitrario
    """

    def __init__(
        self,
        L: int,
        amino_acids: List[str],
        bits_per_pos: int,
        n_qubits: int,
        **kwargs
    ):
        self.L = L
        self.amino_acids = amino_acids
        self.n_aa = len(amino_acids)
        self.bits_per_pos = bits_per_pos
        self.n_qubits = n_qubits
        self.kwargs = kwargs
        self.stats_cache_version = 4

        # 1. Propiedades físico-químicas (Tablas 1 y 2 del paper)
        self._init_properties()

        # 2. Matrices MJ y vecinos estadísticos k=1,3,4
        self._load_matrices_from_files()

        # 3. Calibración Z-score real (Goldstein et al. 1992)
        #    Evalúa n_decoys secuencias aleatorias para estimar μ y σ de cada término
        n_decoys = self.kwargs.get('n_decoys', 10000)
        self.stats = self._load_or_calibrate_z_scores(n_decoys=n_decoys)

    def _default_stats_cache_dir(self) -> str:
        base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".zscore_cache"))
        os.makedirs(base_dir, exist_ok=True)
        return base_dir

    def _stats_cache_payload(self, n_decoys: int) -> Dict[str, Any]:
        return {
            "stats_cache_version": self.stats_cache_version,
            "matrix_hashes": self.matrix_hashes,
            "L": self.L,
            "amino_acids": self.amino_acids,
            "bits_per_pos": self.bits_per_pos,
            "n_qubits": self.n_qubits,
            "n_decoys": int(n_decoys),
            "zscore_seed": int(self.kwargs.get('zscore_seed', 42)),
            "membrane_mode": self.kwargs.get('membrane_mode', 'wheel'),
            "membrane_span": self.kwargs.get('membrane_span'),
            "membrane_positions": self.kwargs.get('membrane_positions'),
            "wheel_phase_deg": float(self.kwargs.get('wheel_phase_deg', 0.0)),
            "wheel_halfwidth_deg": float(self.kwargs.get('wheel_halfwidth_deg', 90.0)),
            "membrane_charge": self.kwargs.get('membrane_charge', 'neg'),
            "max_interaction_dist": int(self.kwargs.get('max_interaction_dist', 1)),
            "electrostatic_cutoff": int(self.kwargs.get('electrostatic_cutoff', 8)),
            "hydro": self.hydro,
            "charges": self.charges,
            "h_alpha": self.h_alpha,
        }

    def _stats_cache_path(self, n_decoys: int) -> str:
        payload = self._stats_cache_payload(n_decoys)
        cache_key = hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()
        cache_dir = self.kwargs.get('stats_cache_dir') or self._default_stats_cache_dir()
        os.makedirs(cache_dir, exist_ok=True)
        return os.path.join(cache_dir, f"zscore_stats_{cache_key}.json")

    def _load_or_calibrate_z_scores(self, n_decoys: int) -> Dict[str, Tuple[float, float]]:
        cache_path = self._stats_cache_path(n_decoys)
        force_recompute = bool(self.kwargs.get('recompute_zscore', False))

        if not force_recompute and os.path.exists(cache_path):
            with open(cache_path, "r", encoding="utf-8") as fh:
                cached = json.load(fh)
            stats = cached.get("stats", {})
            if stats:
                print(f"   ♻️  Reusing cached Z-score calibration: {cache_path}")
                return {
                    key: (float(values["mu"]), float(values["sigma"]))
                    for key, values in stats.items()
                }

        print(f"   📏 Calibrating Z-score statistics with {n_decoys} decoys...")
        stats = self._calibrate_z_scores(n_decoys=n_decoys)
        serializable = {
            key: {"mu": float(mu), "sigma": float(sigma)}
            for key, (mu, sigma) in stats.items()
        }
        payload = {
            "metadata": self._stats_cache_payload(n_decoys),
            "stats": serializable,
        }
        with open(cache_path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2, sort_keys=True)
        print(f"   💾 Saved Z-score calibration cache: {cache_path}")
        return stats

    # ------------------------------------------------------------------
    # Inicialización de propiedades
    # ------------------------------------------------------------------

    def _init_properties(self):
        """
        Pace-Scholtz (hélice) y Fauchère-Pliska (hidrofobicidad).
        Tabla 1 y Tabla 2 del paper.
        """
        # Propensión a hélice (kcal/mol) — menor valor = mayor propensión
        self.h_alpha = {
            'A': 0.00, 'L': 0.21, 'R': 0.21, 'M': 0.24, 'K': 0.26, 'Q': 0.39,
            'E': 0.40, 'I': 0.41, 'W': 0.49, 'S': 0.50, 'Y': 0.53, 'F': 0.54,
            'V': 0.61, 'H': 0.61, 'N': 0.65, 'T': 0.66, 'C': 0.68, 'D': 0.69,
            'G': 1.00, 'P': 3.15
        }
        # Hidrofobicidad Fauchère-Pliska: >0 Apolar, <0 Polar
        self.hydro = {
            'D': -0.77, 'E': -0.64, 'K': -0.99, 'R': -1.01, 'H':  0.13,
            'G':  0.00, 'A':  0.31, 'V':  1.22, 'L':  1.70, 'I':  1.80,
            'P':  0.72, 'M':  1.23, 'F':  1.79, 'W':  2.25, 'Y':  0.96,
            'S': -0.04, 'T':  0.26, 'C':  1.54, 'N': -0.60, 'Q': -0.22
        }
        # Cargas formales
        self.charges = {aa: 0 for aa in self.amino_acids}
        for aa in ['R', 'K']:
            self.charges[aa] = 1
        for aa in ['D', 'E']:
            self.charges[aa] = -1

    def _load_matrices_from_files(self):
        """Map every labelled source matrix into the requested alphabet order."""
        if not self.amino_acids or len(set(self.amino_acids)) != self.n_aa:
            raise ValueError("The amino-acid alphabet must be nonempty and unique")
        self.matrix_hashes = {}
        loaders = {
            "mj_matrix": _load_energy_matrix_file,
            "M1": _load_first_neighbors_matrix_file,
            "M3": _load_third_neighbors_matrix_file,
            "M4": _load_fourth_neighbors_matrix_file,
        }
        for name, loader in loaders.items():
            matrix, symbols = loader()
            aa_to_idx = {aa: i for i, aa in enumerate(symbols)}
            unknown = set(self.amino_acids) - set(symbols)
            if unknown:
                raise ValueError(f"Unknown residues in amino-acid alphabet: {sorted(unknown)}")
            indices = [aa_to_idx[aa] for aa in self.amino_acids]
            setattr(self, name, matrix[np.ix_(indices, indices)].copy())

            # Include the full labelled numerical table in cache identity, not
            # just the selected subset. Explicit little endian gives a stable
            # digest across platforms. MJ is hashed after symmetrization.
            digest = hashlib.sha256(" ".join(symbols).encode("ascii"))
            digest.update(np.asarray(matrix, dtype="<f8").tobytes(order="C"))
            self.matrix_hashes[name] = digest.hexdigest()

    # ------------------------------------------------------------------
    # Entorno (Helical Wheel)
    # ------------------------------------------------------------------

    def _get_env(self, i: int) -> str:
        phi0 = self.kwargs.get('wheel_phase_deg', 0.0)
        # Extraer el ancho de la cara de la membrana (por defecto 90 si no viene)
        half_width = self.kwargs.get('wheel_halfwidth_deg', 90.0)

        # Ángulo del residuo i (100° por residuo en hélice alfa)
        theta = ((i - 1) * 100.0 + phi0) % 360.0

        # Normalizar a rango [-180, 180]
        angle = theta if theta <= 180.0 else theta - 360.0

        return "membrane" if abs(angle) <= half_width else "water"

    # ------------------------------------------------------------------
    # Z-score real  (Sec. 2 del paper: Goldstein 1992)
    # ------------------------------------------------------------------

    def _raw_helix_local(self, seq: List[int]) -> float:
        scale = 1.0 / self.L
        return sum(self.h_alpha.get(self.amino_acids[a], 1.0) for a in seq) * scale

    def _raw_env_pol(self, seq: List[int]) -> float:
        """
        Actua en AMBOS entornos (simetria anfipática):
          Membrana: -h -> premia apolares (h>0), penaliza polares (h<0)
          Agua:     +h -> premia polares  (h<0), penaliza apolares (h>0)
        """
        scale = 1.0 / self.L
        total = 0.0
        for i, a in enumerate(seq, start=1):
            h = self.hydro[self.amino_acids[a]]
            total += -h if self._get_env(i) == "membrane" else +h
        return total * scale

    def _raw_env_chg(self, seq: List[int]) -> float:
        """
        Actua SOLO en agua. Premio UNIDIRECCIONAL: solo favorece la carga
        OPUESTA a sigma, nunca penaliza la carga del mismo signo.
          neg (sigma=-1): K,R,H en agua -> min((-1)(+1), 0) = -1  FAV
                          D,E   en agua -> min((-1)(-1), 0) =  0  NEUTRO
          pos (sigma=+1): D,E   en agua -> min((+1)(-1), 0) = -1  FAV
                          K,R,H en agua -> min((+1)(+1), 0) =  0  NEUTRO
        Residuos neutros (q=0): sin efecto en cualquier caso.
        """
        scale = 1.0 / self.L
        m_charge = self.kwargs.get('membrane_charge', 'neg').lower()
        sigma = {'neg': -1.0, 'pos': 1.0, 'neu': 0.0}.get(m_charge, 0.0)
        total = 0.0
        for i, a in enumerate(seq, start=1):
            if self._get_env(i) != "membrane":   # solo en agua
                raw = sigma * self.charges[self.amino_acids[a]]
                total += min(raw, 0.0)            # solo el componente favorable
        return total * scale

    def _raw_pw_int(self, seq: List[int]) -> float:
        d_max = self.kwargs.get('max_interaction_dist', 1)
        scale = 1.0 / (self.L * d_max)
        total = 0.0
        for i in range(self.L):
            for j in range(i + 1, min(i + d_max + 1, self.L)):
                total += self.mj_matrix[seq[i], seq[j]]
        return total * scale

    def _raw_hydro_moment(self, seq: List[int]) -> float:
        """
        |μ_H|² = (Σ H_α cos(δi) x_{i,α})² + (Σ H_α sin(δi) x_{i,α})²
        Escalado global 2/(L(L-1)).
        """
        scale = 2.0 / (self.L * (self.L - 1)) if self.L > 1 else 1.0
        delta = np.deg2rad(100.0)
        mu_x = sum(self.hydro[self.amino_acids[seq[i]]] * np.cos(delta * (i + 1))
                   for i in range(self.L))
        mu_y = sum(self.hydro[self.amino_acids[seq[i]]] * np.sin(delta * (i + 1))
                   for i in range(self.L))
        return -(mu_x ** 2 + mu_y ** 2) * scale

    def _raw_electrostatic(self, seq: List[int]) -> float:
        r_c = self.kwargs.get('electrostatic_cutoff', 8)
        scale = 2.0 / (self.L * (self.L - 1)) if self.L > 1 else 1.0
        total = 0.0
        for i in range(self.L):
            env_i = self._get_env(i + 1)
            for j in range(i + 1, min(self.L, i + r_c + 1)):
                env_j = self._get_env(j + 1)
                if env_i == "membrane" and env_j == "membrane":
                    eps = 0.1
                elif env_i == "water" and env_j == "water":
                    eps = 1.0
                else:
                    eps = 0.3
                dist = j - i
                qi = self.charges[self.amino_acids[seq[i]]]
                qj = self.charges[self.amino_acids[seq[j]]]
                total += (qi * qj) / ((1 + dist) * eps)
        return total * scale

    def _raw_helix_neigh(self, seq: List[int]) -> float:
        matrices = {1: self.M1, 3: self.M3, 4: self.M4}
        total = 0.0
        for k, matrix in matrices.items():
            if self.L <= k:
                continue
            scale = 1.0 / (self.L - k)
            for i in range(self.L - k):
                total += matrix[seq[i], seq[i + k]] * scale
        return total

    def _calibrate_z_scores(self, n_decoys: int) -> Dict[str, Tuple[float, float]]:
        """
        Estima μ y σ para cada término evaluando n_decoys secuencias aleatorias.
        Implementación fiel al paper (Goldstein et al. 1992, Sec. 2.2).
        """
        rng = np.random.default_rng(self.kwargs.get('zscore_seed', 42))

        accum = {
            'helix_local':   [],
            'env_pol':       [],
            'env_chg':       [],
            'pw_int':        [],
            'hydro_moment':  [],
            'electrostatic': [],
            'helix_neigh':   [],
        }

        for _ in range(n_decoys):
            seq = rng.integers(0, self.n_aa, size=self.L).tolist()
            accum['helix_local'].append(self._raw_helix_local(seq))
            accum['env_pol'].append(self._raw_env_pol(seq))
            accum['env_chg'].append(self._raw_env_chg(seq))
            accum['pw_int'].append(self._raw_pw_int(seq))
            accum['hydro_moment'].append(self._raw_hydro_moment(seq))
            accum['electrostatic'].append(self._raw_electrostatic(seq))
            accum['helix_neigh'].append(self._raw_helix_neigh(seq))

        stats = {}
        for key, values in accum.items():
            arr = np.array(values)
            mu = float(np.mean(arr))
            sigma = float(np.std(arr))
            # Evitar división por cero: si σ≈0 el término es constante → no aporta gradiente
            stats[key] = (mu, sigma if sigma > 1e-12 else 1.0)

        return stats

    def _apply_z_score(self, term_name: str, raw_energy: float) -> float:
        """H' = (H - μ) / σ  (Ec. Z-score del paper)."""
        mu, sigma = self.stats.get(term_name, (0.0, 1.0))
        return (raw_energy - mu) / sigma
