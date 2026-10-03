"""
L-0005: Enstrophy bound for hard spectral Galerkin truncation to |k| ≤ K0.

Let P_{K0} project onto Fourier modes with |k| ≤ K0 (and dealias mask if used).
The Galerkin ODE
  ∂t u = P_{K0} [ -P(u·∇u) + ν Δu ]
keeps support in the shell for all time.

Then L-0001 / L-0003 apply with (K, M) = shell statistics, yielding an explicit
a priori bound. In particular, for t ≤ T,
  Ω(t) ≤ L0001_short(K0)   (often tighter than L0003 for small T).

This PROVES a finite-dimensional bound for the *truncated* dynamics.
It does NOT control the full dealias pseudospectral flow with shell ICs
(C-S-0001), because that flow leaves the shell.

Evidence: N7 (composition of L-0001/L-0004). Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0004_spectral_support import shell_a_priori_bounds
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared, wave_number_grids
from ns_exploration.spectral.integrators import step_etd_rk2, step_rk4, step_semi_implicit_euler
from ns_exploration.spectral.leray import leray_project_hat
from ns_exploration.diagnostics.metrics import compute_diagnostics


def shell_mask(n: int, k_shell: float, dealias: bool = True) -> np.ndarray:
    kr = np.sqrt(k_squared(n))
    mask = kr <= float(k_shell) + 1e-12
    if dealias:
        mask &= dealias_mask(n)
    return mask


def project_to_shell(u_hat: np.ndarray, k_shell: float, dealias: bool = True) -> np.ndarray:
    mask = shell_mask(u_hat.shape[-1], k_shell, dealias=dealias)
    out = u_hat.copy()
    out[:, ~mask] = 0.0
    return leray_project_hat(out)


def step_galerkin_shell(
    u_hat: np.ndarray,
    dt: float,
    nu: float,
    k_shell: float,
    integrator: str = "etd_rk2",
    dealias: bool = True,
) -> np.ndarray:
    """One time step of shell-truncated Galerkin NS."""
    step = {
        "etd_rk2": step_etd_rk2,
        "rk4": step_rk4,
        "semi_implicit": step_semi_implicit_euler,
    }[integrator]
    u = project_to_shell(u_hat, k_shell, dealias=dealias)
    u = step(u, dt, nu, dealias=dealias)
    return project_to_shell(u, k_shell, dealias=dealias)


def evolve_galerkin_shell(
    u_hat: np.ndarray,
    nu: float,
    dt: float,
    t_end: float,
    k_shell: float,
    integrator: str = "etd_rk2",
) -> np.ndarray:
    u = project_to_shell(u_hat, k_shell)
    nsteps = int(round(t_end / dt))
    for _ in range(nsteps):
        u = step_galerkin_shell(u, dt, nu, k_shell, integrator=integrator)
    return u


@dataclass
class LemmaL0005:
    lemma_id: str = "L-0005"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    k_shell: float = 4.0
    n: int = 12
    E0: float = 0.5
    nu: float = 0.1
    t: float = 0.02
    K: float = 0.0
    M: int = 0
    Omega_bound_short: float = 0.0
    Omega_bound_uniform: float = 0.0
    Omega_bound: float = 0.0
    clay_implication: str = (
        "None. Bound for shell-truncated Galerkin ODE only; not full NS on T^3."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0005(
    n: int = 12,
    k_shell: float = 4.0,
    E0: float = 0.5,
    t: float = 0.02,
    nu: float = 0.1,
) -> LemmaL0005:
    b = shell_a_priori_bounds(n, int(k_shell), E0=E0, t=t, nu=nu)
    return LemmaL0005(
        k_shell=k_shell,
        n=n,
        E0=E0,
        nu=nu,
        t=t,
        K=b.K,
        M=b.M,
        Omega_bound_short=b.L0001_short,
        Omega_bound_uniform=b.L0003_uniform,
        Omega_bound=b.best_cap,
        notes=(
            f"For Galerkin truncation |k|≤{k_shell} on N={n}: "
            f"Ω(t)≤{b.best_cap:.6e} for all t∈[0,{t}] (via L-0001∩L-0004). "
            f"K={b.K:.6g}, M={b.M}."
        ),
    )


def save_l0005(lem: LemmaL0005, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(lem.as_dict(), indent=2), encoding="utf-8")


def numerical_check_below_l0005(
    u_hat: np.ndarray,
    k_shell: float,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
) -> dict:
    """N2 check: truncated evolution stays under L-0005 ceiling."""
    n = u_hat.shape[-1]
    lem = lemma_l0005(n=n, k_shell=k_shell, t=t_end, nu=nu)
    uT = evolve_galerkin_shell(u_hat, nu, dt, t_end, k_shell)
    omega = compute_diagnostics(uT, nu).enstrophy
    # verify support
    mask = shell_mask(n, k_shell)
    outside = float(np.sqrt(np.sum(np.abs(uT[:, ~mask]) ** 2)))
    return {
        "Omega_T": omega,
        "L0005_bound": lem.Omega_bound,
        "below_bound": bool(omega <= lem.Omega_bound + 1e-9),
        "energy_outside_shell": outside,
        "lemma": lem.as_dict(),
    }
