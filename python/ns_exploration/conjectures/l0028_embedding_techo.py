"""
L-0028: Embedding techo for L-0027's stretch constant C_†.

L-0027 closes the low-Ω0 slab of C-0007 once
  |⟨ω, curl N⟩| ≤ C Ω^{3/2}  with  C ≤ C_† ≈ 9.562.
The classical Galerkin embedding (L-0001/L-0002)
  ‖∇u‖_∞ ≤ √(3M) √(2Ω)  ⇒  C = 2√2 √(3M)
meets C≤C_† iff M ≤ M_† ≈ 3.81. No dealias subnet of N=24 with a
usable 3D shell has M≤3.81 (the |k|²=1 shell already has 6 modes;
hard truncation |k|²≤2 gives M=18 and C≈20.79 > C_†).

Sobolev/Cauchy variants (Σ|k|^{2−2α} weights, Bernstein–Poincaré) stay
C≳164 on the full mask. Compatible two-point spectrum majorants feeding
the defect N-bound do not push the coupled Ω(T) below M either
(Om0=E0 endpoint ≈41.82 > M).

Hence C_† cannot be reached by L∞ / mode-count embeddings. Closing the
low slab needs a triad-aware stretch bound (cubic-tensor flattening /
shell SOS), not a sharper √M factor. Empirics remain C_emp≪C_† (N2).

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (arithmetic) + N2 (optional stretch probe).
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0002_viscous_galerkin import stretch_constant
from ns_exploration.conjectures.l0023_cubic_dissipation import empirical_stretch_C
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared

def M_star_for_Cdagger(C_dagger: float) -> float:
    """M with 2√2 √(3M) = C_†  ⇒  M = C_†² / 24."""
    return float(C_dagger) ** 2 / 24.0


def modes_upto_radius(n: int, R: int) -> int:
    mask = dealias_mask(n)
    k2 = k_squared(n)
    ctr = Counter(np.rint(k2[mask & (k2 > 0)]).astype(int).tolist())
    return int(sum(m for r, m in ctr.items() if r <= R))


@dataclass
class GalerkinBoundL0028:
    lemma_id: str = "L-0028"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    C_dagger: float = 0.0
    M_star: float = 0.0
    M_full: int = 0
    C_full_L0002: float = 0.0
    M_shell_r1: int = 0
    C_shell_r1: float = 0.0
    M_hard_R2: int = 0
    C_hard_R2: float = 0.0
    embedding_can_meet_Cdagger: bool = False
    C_emp_max: float = 0.0
    emp_below_Cdagger: bool = False
    proved_closes_c0007: bool = False
    clay_implication: str = (
        "None. Embedding techo for stretch C only; C_† still unproved; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0028(
    n: int = 24,
    empirical: bool = True,
    n_random: int = 12,
    ascent_steps: int = 8,
) -> GalerkinBoundL0028:
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    Mstar = M_star_for_Cdagger(Cd)
    M_full = modes_upto_radius(n, 10**9)
    C_full = stretch_constant(M_full)
    M_r1 = modes_upto_radius(n, 1)
    C_r1 = stretch_constant(M_r1)
    M_R2 = modes_upto_radius(n, 2)
    C_R2 = stretch_constant(M_R2)
    C_emp = 0.0
    if empirical:
        emp = empirical_stretch_C(
            n=n, E0=0.5, n_random=n_random, ascent_steps=ascent_steps
        )
        C_emp = float(emp["C_emp_max"])
    can = Mstar >= M_r1 and C_r1 <= Cd + 1e-12  # False in practice
    notes = (
        f"C_†={Cd:.6g} needs M≤M_★={Mstar:.4f} for L-0002 embedding. "
        f"Full mask M={M_full}, C={C_full:.4f}; r=1 shell M={M_r1}, C={C_r1:.4f}; "
        f"|k|²≤2 M={M_R2}, C={C_R2:.4f}. Embedding meets C_†: {can}. "
        f"C_emp≈{C_emp:.6g}. Triad-aware bound required for low-slab closure."
    )
    return GalerkinBoundL0028(
        n=n,
        C_dagger=Cd,
        M_star=Mstar,
        M_full=M_full,
        C_full_L0002=C_full,
        M_shell_r1=M_r1,
        C_shell_r1=C_r1,
        M_hard_R2=M_R2,
        C_hard_R2=C_R2,
        embedding_can_meet_Cdagger=can,
        C_emp_max=C_emp,
        emp_below_Cdagger=C_emp <= Cd + 1e-15,
        proved_closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0028(
    bound: GalerkinBoundL0028,
    path: str | Path = "conjectures/proved_restricted/L-0028.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
