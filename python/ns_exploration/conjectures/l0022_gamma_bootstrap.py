"""
L-0022: Gamma / N_* bootstrap arithmetic for C-0005 (via L-0020).

Conditional (proved): if every dealias field obeys
  ‖N(u)‖_ℓ₂ ≤ γ √(E Ω)
and E≤E0, then under bootstrap Ω≤M on [0,T],
  ‖N‖ ≤ γ √(E0 M).
L-0020 then gives Ω(T)≤M iff γ ≤ γ_crit := N_crit(M)/√(E0 M).

At N=24, E0=0.5, T=0.02, M≈61.237: γ_crit≈1.7469, N_crit≈9.666.

Empirical (N2, VJP ascent on ‖N‖ at fixed energy): adversarial fields reach
  ‖N‖ ≳ 35–40 ≫ N_crit and γ ≳ 10 ≫ γ_crit,
while Ω(T) on those fields (and random/adjoint searches) stays ≲ 41 ≪ M.
Hence the *uniform* N_* / γ hypotheses needed by this bootstrap are FALSE
for all-IC, and L-0020+L-0022 cannot close C-0005 that way.

C-0005 remains open; uniform Duhamel majorants are too crude versus observed
‖N‖. Prior prose “‖N‖≲2” referred only to unstructured random ICs.

FINITE Galerkin / dealias only. Not continuum. Not Clay.
Evidence: N7 (arithmetic) + N2 (adversarial ‖N‖ / Ω(T) probes).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0020_duhamel_h1 import lemma_l0020, omega_duhamel_H1
from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.optimization.adjoint import nonlinear_vjp, real_inner
from ns_exploration.spectral.fourier_conventions import (
    enforce_reality,
    kinetic_energy_from_hat,
)
from ns_exploration.spectral.integrators import step_etd_rk2
from ns_exploration.spectral.leray import leray_project_hat
from ns_exploration.spectral.operators import nonlinear_hat


def gamma_crit(E0: float, M: float, Ncrit: float) -> float:
    return float(Ncrit) / math.sqrt(float(E0) * float(M))


def n_norm(u_hat: np.ndarray) -> float:
    return float(np.sqrt(np.sum(np.abs(nonlinear_hat(u_hat, dealias=True)) ** 2)))


def field_gamma(u_hat: np.ndarray, nu: float = 0.1) -> tuple[float, float, float]:
    e = kinetic_energy_from_hat(u_hat)
    om = compute_diagnostics(u_hat, nu).enstrophy
    nn = n_norm(u_hat)
    if e <= 0.0 or om <= 0.0:
        return nn, om, 0.0
    return nn, om, nn / math.sqrt(e * om)


def _scale_energy(u_hat: np.ndarray, E0: float) -> np.ndarray:
    u_hat = leray_project_hat(enforce_reality(u_hat))
    e = kinetic_energy_from_hat(u_hat)
    if e <= 0:
        raise ValueError("zero energy")
    return u_hat * math.sqrt(E0 / e)


def ascent_N_norm(
    u_hat: np.ndarray,
    E0: float = 0.5,
    steps: int = 40,
    lr: float = 0.08,
) -> tuple[np.ndarray, float]:
    """Riemannian ascent of ‖N‖ on the fixed-energy Leray manifold (N2)."""
    u = _scale_energy(u_hat, E0)
    best_u, best_n = u, n_norm(u)
    for _ in range(steps):
        Nh = nonlinear_hat(u, dealias=True)
        g = nonlinear_vjp(u, Nh, dealias=True)
        g = leray_project_hat(enforce_reality(g))
        uu = real_inner(u, u)
        if uu > 0:
            g = g - (real_inner(g, u) / uu) * u
        gn = math.sqrt(max(real_inner(g, g), 0.0))
        if gn < 1e-16:
            break
        u = _scale_energy(u + (lr / gn) * g, E0)
        nn = n_norm(u)
        if nn > best_n:
            best_n, best_u = nn, u
    return best_u, best_n


def omega_at_T(
    u_hat: np.ndarray,
    nu: float = 0.1,
    dt: float = 1e-3,
    T: float = 0.02,
) -> float:
    u = u_hat.copy()
    for _ in range(int(round(T / dt))):
        u = step_etd_rk2(u, dt, nu, dealias=True)
    return compute_diagnostics(u, nu).enstrophy


def empirical_adversarial_search(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    n_random: int = 12,
    ascent_steps: int = 20,
    seed0: int = 0,
) -> dict:
    """Adversarial ‖N‖ ascent + Ω(T) on those fields (N2)."""
    rng = np.random.default_rng(seed0)
    best = {
        "N_max": 0.0,
        "gamma_max": 0.0,
        "Omega_T_max": 0.0,
        "Omega_T_at_Nmax": 0.0,
    }
    for i in range(n_random):
        k_peak = int(rng.integers(2, 13))
        uh, _ = random_div_free(n, seed=seed0 + 1 + i, k_peak=k_peak, energy_target=E0)
        field, nn = ascent_N_norm(uh, E0=E0, steps=ascent_steps)
        _, om0, g = field_gamma(field, nu=nu)
        omT = omega_at_T(field, nu=nu, T=T)
        if nn > best["N_max"]:
            best["N_max"] = nn
            best["Omega_T_at_Nmax"] = omT
        if g > best["gamma_max"]:
            best["gamma_max"] = g
        if omT > best["Omega_T_max"]:
            best["Omega_T_max"] = omT
        # also raw random Ω(T)
        omT_raw = omega_at_T(uh, nu=nu, T=T)
        if omT_raw > best["Omega_T_max"]:
            best["Omega_T_max"] = omT_raw
    return best


@dataclass
class GalerkinBoundL0022:
    lemma_id: str = "L-0022"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    E0: float = 0.5
    nu: float = 0.1
    T: float = 0.02
    c0005_M: float = 61.23693461895651
    S_T: float = 0.0
    I_sigma: float = 0.0
    Ncrit_c0005: float = 0.0
    gamma_crit: float = 0.0
    N_emp_max: float = 0.0
    gamma_emp_max: float = 0.0
    Omega_T_emp_max: float = 0.0
    Omega_duhamel_at_Nemp: float = 0.0
    uniform_hypothesis_holds_empirically: bool = False
    proved_closes_c0005: bool = False
    clay_implication: str = (
        "None. Bootstrap arithmetic only; uniform γ/N_* false for all-IC; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0022(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    c0005_M: float = 61.23693461895651,
    empirical: bool = True,
    n_random: int = 12,
    ascent_steps: int = 20,
) -> GalerkinBoundL0022:
    b20 = lemma_l0020(n=n, E0=E0, nu=nu, T=T, c0005_M=c0005_M)
    gcrit = gamma_crit(E0, c0005_M, b20.Ncrit_c0005)
    N_emp = 0.0
    g_emp = 0.0
    omT = 0.0
    if empirical:
        emp = empirical_adversarial_search(
            n=n,
            E0=E0,
            nu=nu,
            T=T,
            n_random=n_random,
            ascent_steps=ascent_steps,
        )
        N_emp = float(emp["N_max"])
        g_emp = float(emp["gamma_max"])
        omT = float(emp["Omega_T_max"])
    om_duh = omega_duhamel_H1(b20.S_T, E0, N_emp, b20.I_sigma) if N_emp > 0 else 0.0
    hyp_ok = bool(
        empirical and N_emp <= b20.Ncrit_c0005 + 1e-12 and g_emp <= gcrit + 1e-12
    )
    notes = (
        f"γ_crit={gcrit:.6g}, N_crit={b20.Ncrit_c0005:.6g}. "
        f"Adversarial N2: ‖N‖_max≈{N_emp:.6g}, γ_max≈{g_emp:.6g}, "
        f"Ω(T)_max≈{omT:.6g} (M={c0005_M:.6g}). "
        f"Uniform hypothesis holds empirically: {hyp_ok}. "
        f"Proved C-0005 closure via this bootstrap: False."
    )
    return GalerkinBoundL0022(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        c0005_M=c0005_M,
        S_T=b20.S_T,
        I_sigma=b20.I_sigma,
        Ncrit_c0005=b20.Ncrit_c0005,
        gamma_crit=gcrit,
        N_emp_max=N_emp,
        gamma_emp_max=g_emp,
        Omega_T_emp_max=omT,
        Omega_duhamel_at_Nemp=om_duh,
        uniform_hypothesis_holds_empirically=hyp_ok,
        proved_closes_c0005=False,
        notes=notes,
    )


def save_lemma_l0022(
    bound: GalerkinBoundL0022,
    path: str | Path = "conjectures/proved_restricted/L-0022.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
