"""Time integrators: RK4 and integrating-factor / ETD-RK2. Exploratory (≤ N2)."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np

from ns_exploration.spectral.dealias import apply_dealias
from ns_exploration.spectral.fourier_conventions import k_squared
from ns_exploration.spectral.leray import leray_project_hat
from ns_exploration.spectral.operators import nonlinear_hat


NonlinearFn = Callable[[np.ndarray], np.ndarray]


def rhs_full(u_hat: np.ndarray, nu: float, f_hat: np.ndarray | None, dealias: bool) -> np.ndarray:
    N = nonlinear_hat(u_hat, dealias=dealias)
    from ns_exploration.spectral.operators import viscous_hat

    out = N + viscous_hat(u_hat, nu)
    if f_hat is not None:
        out = out + leray_project_hat(f_hat)
    return out


def step_rk4(
    u_hat: np.ndarray,
    dt: float,
    nu: float,
    f_hat: np.ndarray | None = None,
    dealias: bool = True,
) -> np.ndarray:
    k1 = rhs_full(u_hat, nu, f_hat, dealias)
    k2 = rhs_full(u_hat + 0.5 * dt * k1, nu, f_hat, dealias)
    k3 = rhs_full(u_hat + 0.5 * dt * k2, nu, f_hat, dealias)
    k4 = rhs_full(u_hat + dt * k3, nu, f_hat, dealias)
    out = u_hat + (dt / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
    return leray_project_hat(out)


def step_etd_rk2(
    u_hat: np.ndarray,
    dt: float,
    nu: float,
    f_hat: np.ndarray | None = None,
    dealias: bool = True,
) -> np.ndarray:
    """
    Integrating-factor RK2 (Cox–Matthews style ETD2 approximation for stiff linear term).

    L = -ν|k|². Exact linear propagator φ = exp(L dt).
    """
    k2 = k_squared(u_hat.shape[-1])
    L = -nu * k2
    # Avoid 0/0 at k=0 or ν=0: use series for φ1 = (exp(Ldt)-1)/L
    Ldt = L * dt
    phi = np.exp(Ldt)
    with np.errstate(divide="ignore", invalid="ignore"):
        phi1 = np.where(np.abs(Ldt) < 1e-12, dt * (1 + Ldt / 2), (phi - 1.0) / L)

    def N_only(v: np.ndarray) -> np.ndarray:
        out = nonlinear_hat(v, dealias=dealias)
        if f_hat is not None:
            out = out + leray_project_hat(f_hat)
        return out

    a = N_only(u_hat)
    u_star = phi * u_hat + phi1 * a
    u_star = leray_project_hat(u_star)
    b = N_only(u_star)
    # ETD2: u_{n+1} = φ u + dt*φ1*((a+b)/2) roughly via (φ-1)/L * average
    out = phi * u_hat + phi1 * 0.5 * (a + b)
    return leray_project_hat(out)


def step_semi_implicit_euler(
    u_hat: np.ndarray,
    dt: float,
    nu: float,
    f_hat: np.ndarray | None = None,
    dealias: bool = True,
) -> np.ndarray:
    """Implicit viscosity, explicit nonlinear: û ← (û + dt N) / (1 + dt ν |k|²)."""
    N = nonlinear_hat(u_hat, dealias=dealias)
    if f_hat is not None:
        N = N + leray_project_hat(f_hat)
    denom = 1.0 + dt * nu * k_squared(u_hat.shape[-1])
    out = (u_hat + dt * N) / denom
    return leray_project_hat(out)


INTEGRATORS = {
    "rk4": step_rk4,
    "etd_rk2": step_etd_rk2,
    "semi_implicit": step_semi_implicit_euler,
}
