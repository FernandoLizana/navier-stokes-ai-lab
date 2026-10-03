"""
Polish a candidate IC with adjoint ascent and enforce integrator + resolution gates.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.optimization.adjoint_ascent import adjoint_ascent, enstrophy_integrator
from ns_exploration.optimization.enstrophy_gradient import _project_fixed_energy, _resample_hat
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat
from ns_exploration.spectral.leray import divergence_l2, leray_project_hat


@dataclass
class PolishResult:
    J_start: float
    J_polished_semi: float
    J_etd: float
    J_rk4: float
    J_etd_2n: float
    integrator_rel: float
    resolution_rel: float
    gates_passed: bool
    exceeded_M: bool
    M_target: float
    div_l2: float
    energy: float
    n: int
    evidence_level: str = "N2"
    route: str = "B"

    def as_dict(self) -> dict:
        return asdict(self)


def load_candidate_npz(path: str | Path) -> np.ndarray:
    data = np.load(path)
    return data["u_hat_real"] + 1j * data["u_hat_imag"]


def polish_candidate(
    u_hat: np.ndarray,
    M_target: float,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy_target: float = 0.5,
    n_steps: int = 10,
    step_size: float = 0.06,
    integrator_tol: float = 0.05,
    resolution_tol: float = 0.08,
) -> tuple[np.ndarray, PolishResult]:
    u0 = _project_fixed_energy(leray_project_hat(u_hat), energy_target)
    j0 = enstrophy_integrator(u0, nu, dt, t_end, "etd_rk2")
    u1, asc = adjoint_ascent(
        u0,
        nu=nu,
        dt=dt,
        t_end=t_end,
        energy_target=energy_target,
        n_steps=n_steps,
        step_size=step_size,
        integrator_tol=integrator_tol,
        resolution_tol=resolution_tol,
    )
    j_etd = enstrophy_integrator(u1, nu, dt, t_end, "etd_rk2")
    j_rk4 = enstrophy_integrator(u1, nu, dt, t_end, "rk4")
    integ_rel = abs(j_etd - j_rk4) / (abs(j_etd) + abs(j_rk4) + 1e-30)
    u2 = _project_fixed_energy(_resample_hat(u1, 2 * u1.shape[-1]), energy_target)
    j_2n = enstrophy_integrator(u2, nu, dt, t_end, "etd_rk2")
    res_rel = abs(j_etd - j_2n) / (abs(j_etd) + abs(j_2n) + 1e-30)
    gates = bool(integ_rel <= integrator_tol and res_rel <= resolution_tol)
    # Prefer ascent gates AND our ETD-based checks
    gates = gates and asc.gates_passed
    result = PolishResult(
        J_start=j0,
        J_polished_semi=asc.J_final,
        J_etd=j_etd,
        J_rk4=j_rk4,
        J_etd_2n=j_2n,
        integrator_rel=integ_rel,
        resolution_rel=res_rel,
        gates_passed=gates,
        exceeded_M=bool(j_etd > M_target and gates),
        M_target=M_target,
        div_l2=divergence_l2(u1),
        energy=kinetic_energy_from_hat(u1),
        n=u1.shape[-1],
    )
    return u1, result
