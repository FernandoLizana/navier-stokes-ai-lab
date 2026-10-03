"""
Finite-difference enstrophy gradient at fixed energy (Sprint 2 start).

Evidence ≤ N2. Mandatory resolution-doubling gate before any conjecture.
No claim of maximal enstrophy growth or singularity.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.initial_conditions.generators import random_div_free
from ns_exploration.spectral.fourier_conventions import enforce_reality, kinetic_energy_from_hat
from ns_exploration.spectral.integrators import step_etd_rk2
from ns_exploration.spectral.leray import leray_project_hat


@dataclass
class GradientGateResult:
    n_low: int
    n_high: int
    t_end: float
    nu: float
    energy_target: float
    seed: int
    directional_deriv_low: float
    directional_deriv_high: float
    relative_agreement: float
    gate_passed: bool
    gate_tol: float
    objective_low: float
    objective_high: float
    evidence_level: str = "N2"
    route: str = "B"
    notes: str = "FD directional derivative of enstrophy(T); not adjoint-exact."

    def as_dict(self) -> dict:
        return asdict(self)


def _project_fixed_energy(u_hat: np.ndarray, energy_target: float) -> np.ndarray:
    u_hat = leray_project_hat(enforce_reality(u_hat))
    e = kinetic_energy_from_hat(u_hat)
    if e <= 0:
        raise ValueError("zero energy field")
    return u_hat * np.sqrt(energy_target / e)


def evolve_enstrophy(
    u_hat: np.ndarray,
    nu: float,
    dt: float,
    t_end: float,
) -> float:
    u = u_hat.copy()
    nsteps = int(round(t_end / dt))
    for _ in range(nsteps):
        u = step_etd_rk2(u, dt, nu, dealias=True)
    return compute_diagnostics(u, nu).enstrophy


def _resample_hat(u_hat: np.ndarray, n_new: int) -> np.ndarray:
    """Zero-pad / truncate Fourier coefficients to new grid size."""
    n_old = u_hat.shape[-1]
    out = np.zeros((3, n_new, n_new, n_new), dtype=np.complex128)
    if n_new == n_old:
        return u_hat.copy()

    def copy_octants(src: np.ndarray, dst: np.ndarray) -> None:
        # Copy low modes preserving fft layout (0..pos, then neg)
        no, nn = src.shape[-1], dst.shape[-1]
        h = min(no, nn) // 2
        # positive including 0
        dst[:, :h, :h, :h] = src[:, :h, :h, :h]
        # mixed signs: negative indices
        dst[:, :h, :h, -h:] = src[:, :h, :h, -h:]
        dst[:, :h, -h:, :h] = src[:, :h, -h:, :h]
        dst[:, :h, -h:, -h:] = src[:, :h, -h:, -h:]
        dst[:, -h:, :h, :h] = src[:, -h:, :h, :h]
        dst[:, -h:, :h, -h:] = src[:, -h:, :h, -h:]
        dst[:, -h:, -h:, :h] = src[:, -h:, -h:, :h]
        dst[:, -h:, -h:, -h:] = src[:, -h:, -h:, -h:]

    if n_new > n_old:
        copy_octants(u_hat, out)
    else:
        copy_octants(u_hat, out)
        # when truncating, dst is smaller — copy_octants still works with h=n_new//2
    return leray_project_hat(enforce_reality(out))


def directional_fd_enstrophy_derivative(
    u_hat: np.ndarray,
    direction: np.ndarray,
    nu: float,
    dt: float,
    t_end: float,
    energy_target: float,
    eps: float = 1e-4,
) -> tuple[float, float]:
    """
    Central FD of enstrophy(T) along a divergence-free direction,
    after re-projection to the energy sphere.
    Returns (derivative, objective_at_center).
    """
    d = leray_project_hat(enforce_reality(direction))
    # Tangency to energy sphere: remove radial component in L2(Fourier)
    # <u, d> = Re Σ û* · d̂; project d ← d - (<u,d>/<u,u>) u
    u = _project_fixed_energy(u_hat, energy_target)
    uu = float(np.sum((u.conj() * u).real)) + 1e-30
    ud = float(np.sum((u.conj() * d).real))
    d = d - (ud / uu) * u
    dn = float(np.sqrt(np.sum(np.abs(d) ** 2))) + 1e-30
    d = d / dn

    up = _project_fixed_energy(u + eps * d, energy_target)
    um = _project_fixed_energy(u - eps * d, energy_target)
    jp = evolve_enstrophy(up, nu, dt, t_end)
    jm = evolve_enstrophy(um, nu, dt, t_end)
    j0 = evolve_enstrophy(u, nu, dt, t_end)
    return (jp - jm) / (2 * eps), j0


def resolution_doubling_gate(
    n_low: int = 12,
    seed: int = 0,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    energy_target: float = 0.5,
    gate_tol: float = 0.35,
) -> GradientGateResult:
    """
    Compare directional FD derivative at N and 2N for the same low-mode IC/direction.
    Gate passes if relative agreement ≤ gate_tol.
    """
    n_high = 2 * n_low
    u_low, _ = random_div_free(n_low, seed=seed, energy_target=energy_target)
    # direction: another random field
    d_low, _ = random_div_free(n_low, seed=seed + 17, energy_target=1.0)

    deriv_low, j_low = directional_fd_enstrophy_derivative(
        u_low, d_low, nu, dt, t_end, energy_target
    )

    u_high = _project_fixed_energy(_resample_hat(u_low, n_high), energy_target)
    d_high = _resample_hat(d_low, n_high)
    # Match time step mildly with resolution (keep same dt for Sprint 2 partial)
    deriv_high, j_high = directional_fd_enstrophy_derivative(
        u_high, d_high, nu, dt, t_end, energy_target
    )

    denom = abs(deriv_low) + abs(deriv_high) + 1e-30
    rel = abs(deriv_low - deriv_high) / denom
    return GradientGateResult(
        n_low=n_low,
        n_high=n_high,
        t_end=t_end,
        nu=nu,
        energy_target=energy_target,
        seed=seed,
        directional_deriv_low=deriv_low,
        directional_deriv_high=deriv_high,
        relative_agreement=rel,
        gate_passed=bool(rel <= gate_tol),
        gate_tol=gate_tol,
        objective_low=j_low,
        objective_high=j_high,
    )


def one_gradient_ascent_step(
    u_hat: np.ndarray,
    nu: float,
    dt: float,
    t_end: float,
    energy_target: float,
    step_size: float = 0.05,
    fd_eps: float = 1e-4,
    n_dirs: int = 4,
    seed: int = 0,
) -> tuple[np.ndarray, float, float]:
    """
    Crude multi-directional FD ascent on the energy sphere (low-budget).
    Returns (u_new, J_old, J_new).
    """
    rng = np.random.default_rng(seed)
    u = _project_fixed_energy(u_hat, energy_target)
    j0 = evolve_enstrophy(u, nu, dt, t_end)
    best_dir = None
    best_deriv = -np.inf
    for i in range(n_dirs):
        d = rng.normal(size=u.shape) + 1j * rng.normal(size=u.shape)
        d = leray_project_hat(enforce_reality(d.astype(np.complex128)))
        deriv, _ = directional_fd_enstrophy_derivative(
            u, d, nu, dt, t_end, energy_target, eps=fd_eps
        )
        if deriv > best_deriv:
            best_deriv = deriv
            best_dir = d
    assert best_dir is not None
    # unit tangent direction as in FD helper
    uu = float(np.sum((u.conj() * u).real)) + 1e-30
    ud = float(np.sum((u.conj() * best_dir).real))
    d = best_dir - (ud / uu) * u
    d = d / (float(np.sqrt(np.sum(np.abs(d) ** 2))) + 1e-30)
    u_new = _project_fixed_energy(u + step_size * d, energy_target)
    j1 = evolve_enstrophy(u_new, nu, dt, t_end)
    return u_new, j0, j1
