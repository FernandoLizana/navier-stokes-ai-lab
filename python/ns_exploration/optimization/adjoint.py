"""
Discrete adjoint / VJP for enstrophy(T) on T³ (exploratory ≤ N2).

Forward nonlinear matches operators.nonlinear_hat (convective form + Leray).
Adjoint of DC for C(u)=u·∇u (periodic):
  (DC^* λ)_i = -u·∇λ_i + λ_j ∂_i u_j

Final enstrophy gradient for divergence-free fields:
  ∇_u (½||ω||²) = -Δu  ↔  |k|² û  in Fourier.

Not a Clay claim. Verify against FD directional derivatives.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ns_exploration.spectral.dealias import apply_dealias, dealias_mask
from ns_exploration.spectral.fourier_conventions import (
    fft,
    ifft,
    ifft_vector,
    k_squared,
    wave_number_grids,
)
from ns_exploration.spectral.integrators import step_semi_implicit_euler
from ns_exploration.spectral.leray import leray_project_hat
from ns_exploration.spectral.operators import nonlinear_hat


def real_inner(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.sum((a.conj() * b).real))


def enstrophy_gradient_hat(u_hat: np.ndarray) -> np.ndarray:
    """∇_û (½||ω||²) = |k|² û for div-free fields (Parseval-consistent)."""
    return k_squared(u_hat.shape[-1]) * u_hat


def nonlinear_vjp(u_hat: np.ndarray, lam_hat: np.ndarray, dealias: bool = True) -> np.ndarray:
    """
    VJP of N(u) = -P(u·∇u) at u applied to lam: DN(u)^* lam.
    """
    n = u_hat.shape[-1]
    mask = dealias_mask(n) if dealias else None
    kx, ky, kz = wave_number_grids(n)

    u = ifft_vector(u_hat)
    # Apply P to adjoint seed first (P self-adjoint)
    lam_p = leray_project_hat(lam_hat)
    if dealias:
        lam_p = apply_dealias(lam_p, mask)
    lam = ifft_vector(lam_p)

    # Gradients of u
    du_dx = np.stack([ifft(1j * kx * u_hat[i]) for i in range(3)], axis=0)
    du_dy = np.stack([ifft(1j * ky * u_hat[i]) for i in range(3)], axis=0)
    du_dz = np.stack([ifft(1j * kz * u_hat[i]) for i in range(3)], axis=0)

    # Gradients of lam (use projected spectral adjoint seed)
    dlam_dx = np.stack([ifft(1j * kx * lam_p[i]) for i in range(3)], axis=0)
    dlam_dy = np.stack([ifft(1j * ky * lam_p[i]) for i in range(3)], axis=0)
    dlam_dz = np.stack([ifft(1j * kz * lam_p[i]) for i in range(3)], axis=0)

    # (DC^* λ)_i = -u·∇λ_i + λ_j ∂_i u_j
    adj = np.empty_like(u)
    for i in range(3):
        adv = u[0] * dlam_dx[i] + u[1] * dlam_dy[i] + u[2] * dlam_dz[i]
        # λ_j ∂_i u_j
        if i == 0:
            stretch = lam[0] * du_dx[0] + lam[1] * du_dx[1] + lam[2] * du_dx[2]
        elif i == 1:
            stretch = lam[0] * du_dy[0] + lam[1] * du_dy[1] + lam[2] * du_dy[2]
        else:
            stretch = lam[0] * du_dz[0] + lam[1] * du_dz[1] + lam[2] * du_dz[2]
        adj[i] = -adv + stretch

    adj_hat = np.stack([fft(adj[i]) for i in range(3)], axis=0)
    if dealias:
        adj_hat = apply_dealias(adj_hat, mask)
    # N = -P C ⇒ DN^* = - P ∘ DC^* (with FFT conventions absorbed in phys roundtrip)
    return -leray_project_hat(adj_hat)


def nonlinear_jvp(u_hat: np.ndarray, v_hat: np.ndarray, eps: float = 1e-7, dealias: bool = True) -> np.ndarray:
    """Finite-difference JVP of nonlinear_hat (for cross-checks)."""
    return (nonlinear_hat(u_hat + eps * v_hat, dealias=dealias) - nonlinear_hat(u_hat, dealias=dealias)) / eps


def evolve_store_semi_implicit(
    u_hat: np.ndarray,
    nu: float,
    dt: float,
    t_end: float,
    dealias: bool = True,
) -> list[np.ndarray]:
    """Forward trajectory including u^0 ... u^{N}."""
    u = leray_project_hat(u_hat.copy())
    traj = [u.copy()]
    nsteps = int(round(t_end / dt))
    for _ in range(nsteps):
        u = step_semi_implicit_euler(u, dt, nu, dealias=dealias)
        traj.append(u.copy())
    return traj


def adjoint_gradient_semi_implicit(
    traj: list[np.ndarray],
    nu: float,
    dt: float,
    dealias: bool = True,
) -> np.ndarray:
    """
    Discrete adjoint for u^{n+1} = (u^n + dt N(u^n)) / (1 + dt ν|k|²).

    λ^n = (I + dt DN(u^n)^*) [ λ^{n+1} / (1+dtνk²) ]
    λ^N = |k|² u^N
    Returns λ^0 = ∇_{u^0} enstrophy(T).
    """
    uT = traj[-1]
    lam = enstrophy_gradient_hat(uT)
    k2 = k_squared(uT.shape[-1])
    denom = 1.0 + dt * nu * k2
    for n in range(len(traj) - 2, -1, -1):
        u_n = traj[n]
        lam_over = lam / denom
        vjp = nonlinear_vjp(u_n, lam_over, dealias=dealias)
        lam = lam_over + dt * vjp
        lam = leray_project_hat(lam)
    return lam


def directional_adjoint_derivative(
    u_hat: np.ndarray,
    direction: np.ndarray,
    nu: float,
    dt: float,
    t_end: float,
    dealias: bool = True,
) -> tuple[float, np.ndarray]:
    """⟨∇J, d⟩ using semi-implicit discrete adjoint."""
    traj = evolve_store_semi_implicit(u_hat, nu, dt, t_end, dealias=dealias)
    grad = adjoint_gradient_semi_implicit(traj, nu, dt, dealias=dealias)
    d = leray_project_hat(direction)
    return real_inner(grad, d), grad


def directional_fd_semi_implicit(
    u_hat: np.ndarray,
    direction: np.ndarray,
    nu: float,
    dt: float,
    t_end: float,
    eps: float = 1e-5,
    dealias: bool = True,
) -> float:
    from ns_exploration.diagnostics.metrics import compute_diagnostics

    d = leray_project_hat(direction)
    dn = np.sqrt(real_inner(d, d)) + 1e-30
    d = d / dn
    up = u_hat + eps * d
    um = u_hat - eps * d
    traj_p = evolve_store_semi_implicit(up, nu, dt, t_end, dealias=dealias)
    traj_m = evolve_store_semi_implicit(um, nu, dt, t_end, dealias=dealias)
    jp = compute_diagnostics(traj_p[-1], nu).enstrophy
    jm = compute_diagnostics(traj_m[-1], nu).enstrophy
    return (jp - jm) / (2 * eps)


@dataclass
class AdjointFDAgreement:
    directional_adjoint: float
    directional_fd: float
    relative_error: float
    agreed: bool
    tol: float
    nu: float
    dt: float
    t_end: float
    n: int
    evidence_level: str = "N2"
    route: str = "B"

    def as_dict(self) -> dict:
        return asdict(self)


def verify_adjoint_vs_fd(
    u_hat: np.ndarray,
    direction: np.ndarray,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.01,
    tol: float = 0.15,
) -> AdjointFDAgreement:
    d_adj, _ = directional_adjoint_derivative(u_hat, direction, nu, dt, t_end)
    # FD uses same unit direction scaling as adjoint inner product with raw d —
    # match: normalize d for FD and scale adjoint accordingly
    d = leray_project_hat(direction)
    dn = np.sqrt(real_inner(d, d)) + 1e-30
    d_unit = d / dn
    d_adj_unit, _ = directional_adjoint_derivative(u_hat, d_unit, nu, dt, t_end)
    d_fd = directional_fd_semi_implicit(u_hat, d_unit, nu, dt, t_end)
    rel = abs(d_adj_unit - d_fd) / (abs(d_adj_unit) + abs(d_fd) + 1e-30)
    return AdjointFDAgreement(
        directional_adjoint=d_adj_unit,
        directional_fd=d_fd,
        relative_error=rel,
        agreed=bool(rel <= tol),
        tol=tol,
        nu=nu,
        dt=dt,
        t_end=t_end,
        n=u_hat.shape[-1],
    )
