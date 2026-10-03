"""
L-0023: Cubic-dissipation enstrophy ODE toward C-0005.

Energy is nonincreasing: E(t) ≤ E0. Parseval / Cauchy gives
  Σ_k |k|⁴ ‖û_k‖² ≥ 2 Ω² / E ≥ 2 Ω² / E0.
With N = -P((u·∇)u) and the identity
  dΩ/dt = -ν Σ |k|⁴ ‖û‖² + ⟨ω, curl N⟩,
one obtains
  dΩ/dt ≤ -(2ν/E0) Ω² + ⟨ω, curl N⟩.

If the stretching obeys |⟨ω, curl N⟩| ≤ C Ω^{3/2}, then with z = √Ω,
  z' ≤ -(ν/E0) z³ + (C/2) z².
Starting from the spectral envelope Ω(0) ≤ K²(N) E0, this comparison ODE
yields a finite Ω(T) bound. At N=24, E0=0.5, ν=0.1, T=0.02, K²=147:
  C ≤ C_★ ≈ 2.1539  ⇒  Ω(T) ≤ M_{C-0005} ≈ 61.237.

Empirics (N2): max |stretch|/Ω^{3/2} observed ≪ C_★ (room ~50×).
Proved C ≤ C_★ is still open (L-0002's C ~ √M ≫ C_★).

FINITE Galerkin / dealias only. Not continuum. Not Clay.
Evidence: N7 (ODE arithmetic) + N2 (stretch probes).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.optimization.adjoint import real_inner
from ns_exploration.spectral.fourier_conventions import (
    enforce_reality,
    kinetic_energy_from_hat,
)
from ns_exploration.spectral.leray import leray_project_hat
from ns_exploration.spectral.operators import curl_hat, nonlinear_hat
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def stretch_inner(u_hat: np.ndarray) -> float:
    """⟨ω, curl N⟩ with N = -P((u·∇)u)."""
    w = curl_hat(u_hat)
    n_hat = nonlinear_hat(u_hat, dealias=True)
    return real_inner(w, curl_hat(n_hat))


def stretch_constant_of(u_hat: np.ndarray, nu: float = 0.1) -> tuple[float, float, float]:
    """Return (|stretch|/Ω^{3/2}, Ω, stretch)."""
    om = compute_diagnostics(u_hat, nu).enstrophy
    st = stretch_inner(u_hat)
    if om <= 0:
        return 0.0, om, st
    return abs(st) / (om ** 1.5), om, st


def omega_ode_bound(
    C: float,
    E0: float,
    nu: float,
    T: float,
    Omega0: float,
) -> float:
    """
    Upper bound from z' = -(ν/E0) z³ + (C/2) z², z(0)=√Ω0.
    Returns Ω(T) = z(T)².
    """
    if Omega0 <= 0:
        return 0.0
    z0 = math.sqrt(Omega0)
    alpha = nu / E0
    beta = 0.5 * C

    def rhs(_t: float, z: np.ndarray) -> list[float]:
        zz = float(z[0])
        return [-alpha * zz**3 + beta * zz**2]

    sol = solve_ivp(
        rhs,
        (0.0, float(T)),
        [z0],
        rtol=1e-10,
        atol=1e-12,
        method="RK45",
    )
    if not sol.success:
        raise RuntimeError(f"ODE solve failed: {sol.message}")
    zT = float(sol.y[0, -1])
    return max(zT, 0.0) ** 2


def C_star_for_M(
    E0: float,
    nu: float,
    T: float,
    Omega0: float,
    M: float,
    lo: float = 1e-6,
    hi: float = 20.0,
    iters: int = 60,
) -> float:
    """Largest C with omega_ode_bound(C) ≤ M (binary search)."""
    # Ensure hi is large enough to exceed M
    while omega_ode_bound(hi, E0, nu, T, Omega0) <= M and hi < 1e6:
        hi *= 2.0
    if omega_ode_bound(lo, E0, nu, T, Omega0) > M:
        return 0.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if omega_ode_bound(mid, E0, nu, T, Omega0) <= M:
            lo = mid
        else:
            hi = mid
    return float(lo)


def empirical_stretch_C(
    n: int = 24,
    E0: float = 0.5,
    n_random: int = 40,
    ascent_steps: int = 25,
    seed0: int = 0,
) -> dict:
    """N2 search for C_emp = max |stretch|/Ω^{3/2}."""
    rng = np.random.default_rng(seed0)
    best_C = 0.0
    best_Om = 0.0

    def consider(uh: np.ndarray) -> None:
        nonlocal best_C, best_Om
        c, om, _ = stretch_constant_of(uh)
        if c > best_C:
            best_C, best_Om = c, om

    for i in range(n_random):
        k_peak = int(rng.integers(2, 13))
        uh, _ = random_div_free(n, seed=seed0 + 1 + i, k_peak=k_peak, energy_target=E0)
        consider(uh)
        # random-subspace FD ascent on |stretch|/Ω^{3/2}
        for _ in range(ascent_steps):
            base, _, _ = stretch_constant_of(uh)
            g = np.zeros_like(uh)
            for _dir in range(4):
                v = rng.normal(size=uh.shape) + 1j * rng.normal(size=uh.shape)
                v = leray_project_hat(enforce_reality(v))
                uu = real_inner(uh, uh) + 1e-30
                v = v - (real_inner(v, uh) / uu) * uh
                eps = 1e-4
                up = uh + eps * v
                e = kinetic_energy_from_hat(up)
                if e <= 0:
                    continue
                up = up * math.sqrt(E0 / e)
                d = (stretch_constant_of(up)[0] - base) / eps
                g = g + d * v
            g = leray_project_hat(enforce_reality(g))
            uu = real_inner(uh, uh) + 1e-30
            g = g - (real_inner(g, uh) / uu) * uh
            gn = math.sqrt(max(real_inner(g, g), 0.0))
            if gn < 1e-18:
                break
            uh = uh + (0.06 / gn) * g
            e = kinetic_energy_from_hat(uh)
            uh = leray_project_hat(enforce_reality(uh * math.sqrt(E0 / e)))
            consider(uh)
    return {"C_emp_max": best_C, "Omega_at_Cmax": best_Om}


@dataclass
class GalerkinBoundL0023:
    lemma_id: str = "L-0023"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    E0: float = 0.5
    nu: float = 0.1
    T: float = 0.02
    K2: int = 0
    Omega0_cap: float = 0.0
    c0005_M: float = 61.23693461895651
    C_star: float = 0.0
    Omega_at_Cstar: float = 0.0
    C_emp_max: float = 0.0
    emp_closes: bool = False
    proved_closes_c0005: bool = False
    clay_implication: str = (
        "None. ODE reduction only; stretch constant C≤C_★ not proved; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0023(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    c0005_M: float = 61.23693461895651,
    empirical: bool = True,
    n_random: int = 24,
    ascent_steps: int = 15,
) -> GalerkinBoundL0023:
    K2, _ = full_dealias_exact_stats(n)
    Omega0 = float(K2) * E0
    Cstar = C_star_for_M(E0, nu, T, Omega0, c0005_M)
    om_star = omega_ode_bound(Cstar, E0, nu, T, Omega0)
    C_emp = 0.0
    if empirical:
        emp = empirical_stretch_C(
            n=n, E0=E0, n_random=n_random, ascent_steps=ascent_steps
        )
        C_emp = float(emp["C_emp_max"])
    notes = (
        f"ODE z'=-(ν/E0)z³+(C/2)z² from Ω0≤{Omega0}. "
        f"C_★={Cstar:.6g} ⇒ Ω(T)≤{om_star:.6g} (M={c0005_M:.6g}). "
        f"C_emp≈{C_emp:.6g}; emp_closes={C_emp <= Cstar}. "
        f"Proved C≤C_★: False."
    )
    return GalerkinBoundL0023(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=int(K2),
        Omega0_cap=Omega0,
        c0005_M=c0005_M,
        C_star=Cstar,
        Omega_at_Cstar=om_star,
        C_emp_max=C_emp,
        emp_closes=bool(empirical and C_emp <= Cstar + 1e-12),
        proved_closes_c0005=False,
        notes=notes,
    )


def save_lemma_l0023(
    bound: GalerkinBoundL0023,
    path: str | Path = "conjectures/proved_restricted/L-0023.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
