"""Spectral differential operators and nonlinear term."""

from __future__ import annotations

import numpy as np

from ns_exploration.spectral.dealias import apply_dealias, dealias_mask
from ns_exploration.spectral.fourier_conventions import (
    fft,
    ifft,
    ifft_vector,
    wave_number_grids,
)
from ns_exploration.spectral.leray import leray_project_hat


def curl_hat(u_hat: np.ndarray) -> np.ndarray:
    """ω̂ = ik × û."""
    kx, ky, kz = wave_number_grids(u_hat.shape[-1])
    wx = 1j * (ky * u_hat[2] - kz * u_hat[1])
    wy = 1j * (kz * u_hat[0] - kx * u_hat[2])
    wz = 1j * (kx * u_hat[1] - ky * u_hat[0])
    return np.stack([wx, wy, wz], axis=0)


def grad_scalar_hat(s_hat: np.ndarray) -> np.ndarray:
    kx, ky, kz = wave_number_grids(s_hat.shape[-1])
    return np.stack([1j * kx * s_hat, 1j * ky * s_hat, 1j * kz * s_hat], axis=0)


def nonlinear_hat(
    u_hat: np.ndarray,
    dealias: bool = True,
    mask: np.ndarray | None = None,
) -> np.ndarray:
    """
    Compute -P(u·∇u) in Fourier space via pseudospectral product.

    Uses convective form: (u·∇)u_i = u_j ∂_j u_i.
    """
    n = u_hat.shape[-1]
    if mask is None and dealias:
        mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    u = ifft_vector(u_hat)
    # Derivatives in physical space from spectral
    du_dx = np.stack([ifft(1j * kx * u_hat[i]) for i in range(3)], axis=0)
    du_dy = np.stack([ifft(1j * ky * u_hat[i]) for i in range(3)], axis=0)
    du_dz = np.stack([ifft(1j * kz * u_hat[i]) for i in range(3)], axis=0)
    conv = np.empty_like(u)
    for i in range(3):
        conv[i] = u[0] * du_dx[i] + u[1] * du_dy[i] + u[2] * du_dz[i]
    conv_hat = np.stack([fft(conv[i]) for i in range(3)], axis=0)
    if dealias:
        conv_hat = apply_dealias(conv_hat, mask)
    # RHS of projected NS: -P(conv)
    return -leray_project_hat(conv_hat)


def viscous_hat(u_hat: np.ndarray, nu: float) -> np.ndarray:
    from ns_exploration.spectral.fourier_conventions import k_squared

    return -nu * k_squared(u_hat.shape[-1]) * u_hat
