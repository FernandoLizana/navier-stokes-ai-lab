"""
L-0030: Full-mask triad multiplicity bound on ‖N‖ (toward C_† / C-0007).

Extend the L-0014 counting argument from shell→high to the full nonzero
dealias mask S. With R_★ := max_{k∈S} #{(p,q)∈S² : p+q=k},
  |N̂_k| ≤ Σ_{p+q=k} |q| ‖û_p‖ ‖û_q‖
and Cauchy on each row give
  ‖N‖_ℓ₂ ≤ 2 √(R_★ E Ω).

On N=24 dealias: R_★ = 3148 (vs M=3374 modes). At Ω∼E this beats radial
Young ‖N‖≤2E√ρ_★ (≈78.7 at E=0.5 vs triad ≈56.1). When Ω≫E the triad
bound is weaker than Young; the hybrid
  ‖N‖ ≤ min(2 E √ρ_★, 2 √(R_★ E Ω))
is the natural combination.

Feeding √(2 K² Ω)‖N‖ into the enstrophy comparison on the L-0027 low slab
still yields Ω(T)≫M (hybrid/min-with-defect floors ≳50). The absolute cubic
majorant |stretch|≤(2Ω)^{3/2}√R_★ would give C≤2√2√R_★≈158.7, still ≫C_†≈9.56
and only marginally below L-0029's C_F≈164.3.

So triad *counting* improves ‖N‖ near equipartition but does not reach C_†.
Weighted triad / SOS flattenings remain the next target.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (counting + Cauchy) + N2 (ODE probes).
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ns_exploration.conjectures.l0020_duhamel_h1 import full_mask_rho_star
from ns_exploration.conjectures.l0024_spectral_defect import (
    N_defect_bound,
    dealias_shell_stats,
)
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0029_frobenius_stretch import stretch_constant_frobenius
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def dealias_nonzero_modes(n: int) -> list[tuple[int, int, int]]:
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    return [
        (int(a), int(b), int(c))
        for a, b, c in zip(kx[mask], ky[mask], kz[mask])
        if not (a == 0 and b == 0 and c == 0)
    ]


def full_mask_triad_Rstar(n: int = 24) -> tuple[int, int]:
    """Return (R_★, n_triads) on the nonzero dealias mask."""
    modes = dealias_nonzero_modes(n)
    mset = set(modes)
    ctr: Counter[tuple[int, int, int]] = Counter()
    for p in modes:
        px, py, pz = p
        for qx, qy, qz in modes:
            s = (px + qx, py + qy, pz + qz)
            if s in mset:
                ctr[s] += 1
    if not ctr:
        return 0, 0
    return int(max(ctr.values())), int(sum(ctr.values()))


def N_triad_bound(E: float, Omega: float, R_star: int) -> float:
    """‖N‖ ≤ 2 √(R_★ E Ω)."""
    E = max(float(E), 0.0)
    Omega = max(float(Omega), 0.0)
    return 2.0 * math.sqrt(float(R_star) * E * Omega)


def N_young_bound(E: float, rho_star: float) -> float:
    return 2.0 * float(E) * math.sqrt(float(rho_star))


def N_hybrid_bound(E: float, Omega: float, R_star: int, rho_star: float) -> float:
    return min(N_young_bound(E, rho_star), N_triad_bound(E, Omega, R_star))


def cubic_C_from_Rstar(R_star: int) -> float:
    """C ≤ 2√2 √R_★ under the absolute a-cubic majorant Σ a_p a_q a_s."""
    return 2.0 * math.sqrt(2.0) * math.sqrt(float(R_star))


def low_slab_ode_worst(
    N_fn,
    E0: float,
    nu: float,
    T: float,
    K2: int,
    Omega_hi: float,
    n_grid: int = 15,
) -> float:
    def ode(Om0: float) -> float:
        def rhs(_t: float, y: np.ndarray) -> list[float]:
            E = max(float(y[0]), 1e-15)
            Om = min(max(float(y[1]), 1e-15), float(K2) * E)
            Nbd = float(N_fn(E, Om))
            return [-2.0 * nu * Om, -(2.0 * nu / E) * Om * Om + math.sqrt(2.0 * K2 * Om) * Nbd]

        sol = solve_ivp(
            rhs,
            (0.0, float(T)),
            [float(E0), float(Om0)],
            rtol=1e-8,
            atol=1e-11,
            max_step=1e-4,
        )
        return float(sol.y[1, -1])

    xs = np.linspace(E0, Omega_hi, n_grid)
    return max(ode(float(o)) for o in xs)


@dataclass
class GalerkinBoundL0030:
    lemma_id: str = "L-0030"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    M_modes: int = 0
    R_star: int = 0
    n_triads: int = 0
    rho_star: float = 0.0
    C_dagger: float = 0.0
    C_from_Rstar: float = 0.0
    C_F_L0029: float = 0.0
    N_triad_at_equipartition: float = 0.0
    N_young_at_E0: float = 0.0
    low_slab_hybrid_worst: float = 0.0
    low_slab_min_defect_worst: float = 0.0
    c0007_M: float = 0.0
    closes_c0007: bool = False
    clay_implication: str = (
        "None. Full-mask triad counting for ‖N‖ only; C_† not reached; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0030(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    n_grid: int = 15,
) -> GalerkinBoundL0030:
    b7 = lemma_l0027(n=n, empirical=False)
    R_star, n_tri = full_mask_triad_Rstar(n)
    rho, _, _ = full_mask_rho_star(n)
    K2, _ = full_dealias_exact_stats(n)
    M = len(dealias_nonzero_modes(n))
    stats = dealias_shell_stats(n)
    C_R = cubic_C_from_Rstar(R_star)
    C_F = stretch_constant_frobenius(M)
    N_eq = N_triad_bound(E0, E0, R_star)
    N_y = N_young_bound(E0, rho)

    def hybrid(E: float, Om: float) -> float:
        return N_hybrid_bound(E, Om, R_star, rho)

    def min_def(E: float, Om: float) -> float:
        return min(hybrid(E, Om), N_defect_bound(E, Om, stats))

    hy = low_slab_ode_worst(hybrid, E0, nu, T, int(K2), b7.Omega_star, n_grid=n_grid)
    md = low_slab_ode_worst(min_def, E0, nu, T, int(K2), b7.Omega_star, n_grid=n_grid)
    notes = (
        f"R_★={R_star}, n_triads={n_tri}, M={M}. "
        f"‖N‖≤2√(R_★EΩ): at Ω=E={E0} gives {N_eq:.4g} vs Young {N_y:.4g}. "
        f"Cubic C≤2√2√R_★={C_R:.4g} (C_F={C_F:.4g}, C_†={b7.C_dagger:.4g}). "
        f"Low-slab hybrid Ω(T)≤{hy:.4g}, min(hybrid,defect)≤{md:.4g} "
        f"(M={b7.c0007_M:.4g}); closes C-0007: False."
    )
    return GalerkinBoundL0030(
        n=n,
        M_modes=M,
        R_star=R_star,
        n_triads=n_tri,
        rho_star=float(rho),
        C_dagger=b7.C_dagger,
        C_from_Rstar=C_R,
        C_F_L0029=C_F,
        N_triad_at_equipartition=N_eq,
        N_young_at_E0=N_y,
        low_slab_hybrid_worst=hy,
        low_slab_min_defect_worst=md,
        c0007_M=b7.c0007_M,
        closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0030(
    bound: GalerkinBoundL0030,
    path: str | Path = "conjectures/proved_restricted/L-0030.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
