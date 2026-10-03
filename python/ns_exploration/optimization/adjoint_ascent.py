"""
Adjoint gradient ascent of enstrophy(T) on the fixed-energy sphere.

Uses semi-implicit discrete adjoint (optimization/adjoint.py).
Mandatory post-checks: N↔2N objective agreement, RK4↔ETD (exploratory ≤ N2).
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.optimization.adjoint import (
    adjoint_gradient_semi_implicit,
    evolve_store_semi_implicit,
    real_inner,
)
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy, _resample_hat
from ns_exploration.spectral.fourier_conventions import enforce_reality, kinetic_energy_from_hat
from ns_exploration.spectral.integrators import step_etd_rk2, step_rk4, step_semi_implicit_euler
from ns_exploration.spectral.leray import divergence_l2, leray_project_hat


def _tangent_project(u: np.ndarray, g: np.ndarray) -> np.ndarray:
    """Project g onto tangent of energy sphere at u (and Leray)."""
    g = leray_project_hat(enforce_reality(g))
    uu = real_inner(u, u) + 1e-30
    ug = real_inner(u, g)
    g = g - (ug / uu) * u
    return leray_project_hat(g)


def enstrophy_semi_implicit(u_hat: np.ndarray, nu: float, dt: float, t_end: float) -> float:
    traj = evolve_store_semi_implicit(u_hat, nu, dt, t_end)
    return compute_diagnostics(traj[-1], nu).enstrophy


def enstrophy_integrator(u_hat: np.ndarray, nu: float, dt: float, t_end: float, name: str) -> float:
    step = {"etd_rk2": step_etd_rk2, "rk4": step_rk4, "semi_implicit": step_semi_implicit_euler}[name]
    u = u_hat.copy()
    for _ in range(int(round(t_end / dt))):
        u = step(u, dt, nu, dealias=True)
    return compute_diagnostics(u, nu).enstrophy


@dataclass
class AscentResult:
    n: int
    steps: int
    J_history: list[float]
    J_final: float
    J_etd: float
    J_rk4: float
    J_semi_2n: float
    integrator_rel: float
    resolution_rel: float
    gates_passed: bool
    div_l2: float
    energy: float
    evidence_level: str = "N2"
    route: str = "B"

    def as_dict(self) -> dict:
        return asdict(self)


def adjoint_ascent(
    u0: np.ndarray,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy_target: float = 0.5,
    n_steps: int = 8,
    step_size: float = 0.05,
    integrator_tol: float = 0.05,
    resolution_tol: float = 0.08,
) -> tuple[np.ndarray, AscentResult]:
    """
    Gradient ascent on energy sphere using discrete adjoint direction.
    Line search: backtrack if objective does not increase.
    """
    u = _project_fixed_energy(u0, energy_target)
    history: list[float] = []
    j = enstrophy_semi_implicit(u, nu, dt, t_end)
    history.append(j)

    for _ in range(n_steps):
        traj = evolve_store_semi_implicit(u, nu, dt, t_end)
        grad = adjoint_gradient_semi_implicit(traj, nu, dt)
        direction = _tangent_project(u, grad)
        gn = np.sqrt(real_inner(direction, direction)) + 1e-30
        direction = direction / gn

        ss = step_size
        improved = False
        for _bt in range(8):
            u_trial = _project_fixed_energy(u + ss * direction, energy_target)
            j_trial = enstrophy_semi_implicit(u_trial, nu, dt, t_end)
            if j_trial > j + 1e-10:
                u = u_trial
                j = j_trial
                improved = True
                break
            ss *= 0.5
        history.append(j)
        if not improved:
            break

    # Gates: integrator cross-check at N; resolution via semi-implicit at 2N
    j_etd = enstrophy_integrator(u, nu, dt, t_end, "etd_rk2")
    j_rk4 = enstrophy_integrator(u, nu, dt, t_end, "rk4")
    integ_rel = abs(j_etd - j_rk4) / (abs(j_etd) + abs(j_rk4) + 1e-30)

    u2 = _project_fixed_energy(_resample_hat(u, 2 * u.shape[-1]), energy_target)
    j_2n = enstrophy_semi_implicit(u2, nu, dt, t_end)
    j_n_semi = enstrophy_semi_implicit(u, nu, dt, t_end)
    res_rel = abs(j_n_semi - j_2n) / (abs(j_n_semi) + abs(j_2n) + 1e-30)

    gates_ok = bool(integ_rel <= integrator_tol and res_rel <= resolution_tol)
    result = AscentResult(
        n=u.shape[-1],
        steps=len(history) - 1,
        J_history=history,
        J_final=j,
        J_etd=j_etd,
        J_rk4=j_rk4,
        J_semi_2n=j_2n,
        integrator_rel=integ_rel,
        resolution_rel=res_rel,
        gates_passed=gates_ok,
        div_l2=divergence_l2(u),
        energy=kinetic_energy_from_hat(u),
    )
    return u, result
