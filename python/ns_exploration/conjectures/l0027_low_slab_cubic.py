"""
L-0027: Low-Ω0 cubic conditional toward C-0007.

L-0026 splits C-0007 at Ω★≈18.87: the high slab Ω0≥Ω★ is already ≤M via
L-0024. Remains the low slab Ω0∈[E0,Ω★).

On that slab, the L-0023 comparison
  z' = -(ν/E0) z³ + (C/2) z²,  z=√Ω,
closes Ω(T)≤M≈41.284 uniformly in Ω0≤Ω★ as soon as
  C ≤ C_† ≈ 9.562
(binary search; worst at Ω0=Ω★). Empirically C_emp≪C_† (N2). Proved
C≤C_† is still open (L-0002 gives C∼√M_modes≃285≫C_†).

Note: the same cubic ODE cannot close *all*-IC C-0007 even at C=0
(majorant from Ω0=K²E0 exceeds M); the high slab must use L-0024/L-0026.

Together: (L-0026 high slab) + (C≤C_† on the dynamics) ⇒ C-0007.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (ODE arithmetic) + N2 (optional stretch probe).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0023_cubic_dissipation import (
    empirical_stretch_C,
    omega_ode_bound,
)
from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def C_star_low_slab(
    E0: float,
    nu: float,
    T: float,
    M: float,
    Omega_star: float,
    n_grid: int = 25,
    lo: float = 1e-6,
    hi: float = 40.0,
    iters: int = 50,
) -> tuple[float, float]:
    """Largest C with max_{Ω0∈[E0,Ω★]} omega_ode_bound(C,Ω0) ≤ M."""

    def worst(C: float) -> float:
        xs = np.linspace(E0, Omega_star, n_grid)
        return max(omega_ode_bound(C, E0, nu, T, float(o)) for o in xs)

    while worst(hi) <= M and hi < 1e5:
        hi *= 2.0
    if worst(lo) > M:
        return 0.0, worst(lo)
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if worst(mid) <= M:
            lo = mid
        else:
            hi = mid
    return float(lo), float(worst(lo))


def all_ic_cubic_floor(
    E0: float,
    nu: float,
    T: float,
    K2: int,
    n_grid: int = 41,
) -> float:
    """Worst Ω(T) of L-0023 ODE at C=0 over Ω0∈[E0,K²E0]."""
    xs = np.linspace(E0, float(K2) * E0, n_grid)
    return max(omega_ode_bound(0.0, E0, nu, T, float(o)) for o in xs)


@dataclass
class GalerkinBoundL0027:
    lemma_id: str = "L-0027"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    E0: float = 0.5
    nu: float = 0.1
    T: float = 0.02
    c0007_M: float = 0.0
    Omega_star: float = 0.0
    C_dagger: float = 0.0
    Omega_T_at_Cdagger: float = 0.0
    C_emp_max: float = 0.0
    emp_closes_low_slab: bool = False
    all_ic_cubic_C0_floor: float = 0.0
    all_ic_cubic_closes_c0007: bool = False
    proved_closes_c0007: bool = False
    clay_implication: str = (
        "None. Low-slab cubic conditional only; C≤C_† not proved; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0027(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    empirical: bool = True,
    n_random: int = 16,
    ascent_steps: int = 10,
    n_grid: int = 25,
) -> GalerkinBoundL0027:
    b26 = lemma_l0026(n=n, E0=E0, nu=nu, T=T, n_grid=61)
    M = b26.c0007_M
    om_star = b26.Omega_star
    Cd, omT = C_star_low_slab(E0, nu, T, M, om_star, n_grid=n_grid)
    K2, _ = full_dealias_exact_stats(n)
    floor0 = all_ic_cubic_floor(E0, nu, T, K2, n_grid=41)
    C_emp = 0.0
    if empirical:
        emp = empirical_stretch_C(
            n=n, E0=E0, n_random=n_random, ascent_steps=ascent_steps
        )
        C_emp = float(emp["C_emp_max"])
    notes = (
        f"Low slab Ω0≤Ω★={om_star:.6f}: C_†={Cd:.6g} ⇒ Ω(T)≤{omT:.6g} (M={M:.6g}). "
        f"C_emp≈{C_emp:.6g}; emp_closes={C_emp <= Cd}. "
        f"All-IC cubic at C=0 already has floor={floor0:.6g}>M "
        f"(high slab needs L-0024/L-0026). Proved C≤C_†: False."
    )
    return GalerkinBoundL0027(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        c0007_M=M,
        Omega_star=om_star,
        C_dagger=Cd,
        Omega_T_at_Cdagger=omT,
        C_emp_max=C_emp,
        emp_closes_low_slab=C_emp <= Cd + 1e-15,
        all_ic_cubic_C0_floor=floor0,
        all_ic_cubic_closes_c0007=floor0 <= M + 1e-9,
        proved_closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0027(
    bound: GalerkinBoundL0027,
    path: str | Path = "conjectures/proved_restricted/L-0027.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
