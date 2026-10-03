"""
Quantitative residuals for known fields (exploratory ≤ N3).

Beltrami (ABC): P((u·∇)u) = 0 exactly in the continuum.
Manufactured linear mode: exact viscous decay with zero nonlinear term.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ns_exploration.initial_conditions.generators import abc_flow
from ns_exploration.spectral.fourier_conventions import kinetic_energy_from_hat, wave_number_grids
from ns_exploration.spectral.integrators import rhs_full, step_etd_rk2
from ns_exploration.spectral.leray import divergence_l2, leray_project_hat
from ns_exploration.spectral.operators import nonlinear_hat, viscous_hat


@dataclass
class ResidualReport:
    name: str
    n: int
    dealias: bool
    residual_l2: float
    residual_rel: float
    divergence_l2: float
    energy: float
    notes: str
    evidence_level: str = "N2"
    route: str = "B"

    def as_dict(self) -> dict:
        return asdict(self)


def _l2_hat(field_hat: np.ndarray) -> float:
    return float(np.sqrt(np.sum(np.abs(field_hat) ** 2)))


def beltrami_nonlinear_residual(n: int, dealias: bool = True) -> ResidualReport:
    """
    For ABC Beltrami field, continuum P((u·∇)u)=0, so nonlinear_hat ≈ 0.
    Reports ||N(u)||_2 and relative to ||u||_2.
    """
    u_hat, _ = abc_flow(n)
    N = nonlinear_hat(u_hat, dealias=dealias)
    # nonlinear_hat returns -P(conv); residual of Beltrami identity is ||P(conv)|| = ||N||
    res = _l2_hat(N)
    scale = _l2_hat(u_hat) + 1e-30
    return ResidualReport(
        name="abc_beltrami_nonlinear",
        n=n,
        dealias=dealias,
        residual_l2=res,
        residual_rel=res / scale,
        divergence_l2=divergence_l2(u_hat),
        energy=kinetic_energy_from_hat(u_hat),
        notes="Continuum expects residual=0; nonzero = dealiasing/FFT roundoff.",
    )


def manufactured_linear_mode_residual(
    n: int,
    nu: float = 0.1,
    dt: float = 1e-3,
    steps: int = 20,
) -> ResidualReport:
    """
    Single div-free Fourier mode: exact solution û(t)=û0 e^{-ν|k|² t}.
    Residual = ||u_num(T) - u_exact(T)|| / ||u_exact(T)|| after ETD steps.
    """
    uh = np.zeros((3, n, n, n), dtype=np.complex128)
    uh[1, 1, 0, 0] = 1e-3
    uh[1, -1, 0, 0] = 1e-3
    uh = leray_project_hat(uh)
    u0 = uh.copy()
    for _ in range(steps):
        uh = step_etd_rk2(uh, dt, nu, dealias=True)
    t = dt * steps
    kx, ky, kz = wave_number_grids(n)
    decay = np.exp(-nu * (kx * kx + ky * ky + kz * kz) * t)
    exact = u0 * decay
    err = _l2_hat(uh - exact)
    scale = _l2_hat(exact) + 1e-30
    return ResidualReport(
        name="manufactured_linear_mode",
        n=n,
        dealias=True,
        residual_l2=err,
        residual_rel=err / scale,
        divergence_l2=divergence_l2(uh),
        energy=kinetic_energy_from_hat(uh),
        notes="Exact viscous decay of one mode; nonlinear negligible at tiny amplitude.",
    )


def rhs_balance_residual(u_hat: np.ndarray, nu: float, dealias: bool = True) -> float:
    """||∂t u_rhs|| diagnostic magnitude (not a PDE residual against data)."""
    return _l2_hat(rhs_full(u_hat, nu, None, dealias))


def compare_residuals_two_resolutions(
    ns: tuple[int, int] = (16, 32),
) -> list[ResidualReport]:
    reports: list[ResidualReport] = []
    for n in ns:
        reports.append(beltrami_nonlinear_residual(n, dealias=True))
        reports.append(beltrami_nonlinear_residual(n, dealias=False))
        reports.append(manufactured_linear_mode_residual(n))
    return reports
