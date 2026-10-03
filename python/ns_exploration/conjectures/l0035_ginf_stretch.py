"""
L-0035: Shell two-point ‖∇u‖_{∞,F} stretch bound (C_† / C-0007 techo).

Pointwise Frobenius:
  ‖∇u(x)‖_F ≤ Σ_k |k| ‖û_k‖ ≤ Σ_r √(r m_r) √(2 e_r) =: g_∞(e).
The map e ↦ g_∞ is concave, hence maximized under Σ e_r=E, Σ r e_r=Ω at
extreme (two-point / singleton) spectra. Then
  |stretch| ≤ 2 Ω g_∞  ⇒  C_eff = 2 g_∞ / √Ω
(compare L-0029 absolute Frobenius, refined by shell two-point).

Consequences on N=24 dealias (E0=0.5, ν=0.1, T=0.02):
- At Ω=E=0.5 (r=1 singleton): C_eff ≈ 6.93 ≤ C_† ≈ 9.56.
- On the low slab Ω∈[E0,Ω★], max C_eff ≳ 36 ≫ C_†; the set where
  C_eff≤C_† is essentially only the equipartition corner Ω≈E.
- Pure g_∞ comparison ODE on the low slab floors ≳80 > M.
- Hybrid production min(g_∞, defect, shell-N, Young) still floors ≳55 > M.

So the shell-refined ‖∇u‖_∞ bound meets C_† only at equipartition and does
not close the low slab, alone or hybridized with prior majorants.
Still missing: signed / SOS triad cancellation.

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (shell estimates + two-point extrema) + N2 (ODE probes).
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
from ns_exploration.conjectures.l0034_shell_N_bound import N_shell_bound, shell_multiplicities
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def shell_radii_masses(n: int = 24) -> tuple[list[int], list[int]]:
    mask = dealias_mask(n)
    k2 = k_squared(n)
    ctr = Counter(np.rint(k2[mask & (k2 > 0)]).astype(int).tolist())
    radii = sorted(ctr.keys())
    return radii, [ctr[r] for r in radii]


def max_ginf_two_point(E: float, Om: float, radii: list[int], masses: list[int]) -> float:
    """Max Σ_r √(r m_r) √(2 e_r) over two-point spectra with given (E,Ω)."""
    best = 0.0
    E = float(E)
    Om = float(Om)
    for i, r1 in enumerate(radii):
        m1 = masses[i]
        for j in range(i, len(radii)):
            r2 = radii[j]
            m2 = masses[j]
            if r1 == r2:
                if abs(Om - r1 * E) <= 1e-9 * max(E, 1.0):
                    best = max(best, math.sqrt(r1 * m1) * math.sqrt(2.0 * E))
                continue
            if Om < r1 * E - 1e-12 or Om > r2 * E + 1e-12:
                continue
            e2 = (Om - r1 * E) / (r2 - r1)
            e1 = E - e2
            if e1 < 0.0 or e2 < 0.0:
                continue
            g = math.sqrt(r1 * m1) * math.sqrt(2.0 * e1) + math.sqrt(r2 * m2) * math.sqrt(
                2.0 * e2
            )
            best = max(best, g)
    return float(best)


def C_eff_ginf(E: float, Om: float, radii: list[int], masses: list[int]) -> float:
    Om = max(float(Om), 1e-15)
    g = max_ginf_two_point(E, Om, radii, masses)
    return 2.0 * g / math.sqrt(Om)


def stretch_ginf(E: float, Om: float, radii: list[int], masses: list[int]) -> float:
    return 2.0 * float(Om) * max_ginf_two_point(E, Om, radii, masses)


def production_min(
    E: float,
    Om: float,
    radii: list[int],
    masses: list[int],
    stats: dict,
    rho_star: float,
    K2: int,
) -> float:
    p_ginf = stretch_ginf(E, Om, radii, masses)
    Nd = N_defect_bound(E, Om, stats)
    p_def = math.sqrt(2.0 * K2 * Om) * Nd
    Ns = N_shell_bound(E, Om, radii, masses, rho_star)
    p_sh = math.sqrt(2.0 * K2 * Om) * Ns
    Ny = 2.0 * E * math.sqrt(rho_star)
    p_y = math.sqrt(2.0 * K2 * Om) * Ny
    return min(p_ginf, p_def, p_sh, p_y)


def _ode_worst(
    E0: float,
    nu: float,
    T: float,
    K2: int,
    Omega_hi: float,
    rhs_prod,
    n_grid: int = 9,
) -> float:
    def ode(Om0: float) -> float:
        def rhs(_t: float, y: np.ndarray) -> list[float]:
            E = max(float(y[0]), 1e-15)
            Om = min(max(float(y[1]), 1e-15), float(K2) * E)
            return [-2.0 * nu * Om, -(2.0 * nu / E) * Om * Om + rhs_prod(E, Om)]

        sol = solve_ivp(
            rhs, (0.0, float(T)), [float(E0), float(Om0)], rtol=1e-7, atol=1e-10, max_step=5e-4
        )
        return float(sol.y[1, -1])

    xs = np.linspace(E0, Omega_hi, n_grid)
    return max(ode(float(o)) for o in xs)


def C_max_on_slab(
    E0: float,
    Omega_hi: float,
    radii: list[int],
    masses: list[int],
    n_grid: int = 41,
) -> tuple[float, float]:
    best_C = 0.0
    best_Om = float(E0)
    for Om in np.linspace(E0, Omega_hi, n_grid):
        C = C_eff_ginf(E0, float(Om), radii, masses)
        if C > best_C:
            best_C = C
            best_Om = float(Om)
    return best_C, best_Om


def Omega_c_for_Cdagger(
    E0: float,
    Omega_hi: float,
    C_dagger: float,
    radii: list[int],
    masses: list[int],
    n_grid: int = 401,
) -> float:
    """Largest Ω in [E0,Ω★] with C_eff(E0,Ω) ≤ C_† (scan; may equal E0)."""
    last = float(E0)
    for Om in np.linspace(E0, Omega_hi, n_grid):
        if C_eff_ginf(E0, float(Om), radii, masses) <= C_dagger + 1e-12:
            last = float(Om)
        else:
            break
    return last


@dataclass
class GalerkinBoundL0035:
    lemma_id: str = "L-0035"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    C_at_equipartition: float = 0.0
    C_max_low_slab: float = 0.0
    Omega_c: float = 0.0
    ginf_ode_worst: float = 0.0
    minprod_ode_worst: float = 0.0
    closes_c0007: bool = False
    clay_implication: str = (
        "None. Shell two-point ‖∇u‖_∞ stretch without triad cancellation; "
        "C_† only at equipartition; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0035(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    n_grid: int = 9,
) -> GalerkinBoundL0035:
    b7 = lemma_l0027(n=n, empirical=False)
    radii, masses = shell_radii_masses(n)
    # reuse L-0034 shell list for N_shell_bound compatibility
    radii34, masses34 = shell_multiplicities(n)
    assert radii == radii34 and masses == masses34
    rho, _, _ = full_mask_rho_star(n)
    K2, _ = full_dealias_exact_stats(n)
    stats = dealias_shell_stats(n)
    C_eq = C_eff_ginf(E0, E0, radii, masses)
    C_max, _ = C_max_on_slab(E0, b7.Omega_star, radii, masses, n_grid=41)
    Om_c = Omega_c_for_Cdagger(E0, b7.Omega_star, b7.C_dagger, radii, masses)

    def prod_ginf(E: float, Om: float) -> float:
        return stretch_ginf(E, Om, radii, masses)

    def prod_min(E: float, Om: float) -> float:
        return production_min(E, Om, radii, masses, stats, rho, int(K2))

    g_ode = _ode_worst(E0, nu, T, int(K2), b7.Omega_star, prod_ginf, n_grid=n_grid)
    m_ode = _ode_worst(E0, nu, T, int(K2), b7.Omega_star, prod_min, n_grid=n_grid)
    notes = (
        f"Shell g_∞ two-point: C(E)={C_eq:.4g}≤C_†={b7.C_dagger:.4g}; "
        f"C_max(low)={C_max:.4g}; Ω_c≈{Om_c:.4g}. "
        f"ginf ODE worst={g_ode:.4g}, min-prod ODE worst={m_ode:.4g} "
        f"(M={b7.c0007_M:.4g}). Closes C-0007: False."
    )
    return GalerkinBoundL0035(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=b7.C_dagger,
        C_at_equipartition=C_eq,
        C_max_low_slab=C_max,
        Omega_c=Om_c,
        ginf_ode_worst=g_ode,
        minprod_ode_worst=m_ode,
        closes_c0007=False,
        notes=notes,
    )


def save_lemma_l0035(
    bound: GalerkinBoundL0035,
    path: str | Path = "conjectures/proved_restricted/L-0035.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
