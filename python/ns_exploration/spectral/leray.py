"""Leray projector in Fourier space. Route B."""

from __future__ import annotations

import numpy as np

from ns_exploration.spectral.fourier_conventions import wave_number_grids


def leray_project_hat(u_hat: np.ndarray) -> np.ndarray:
    """
    Apply P_k = I - kkᵀ/|k|² modewise. Zero mode unchanged.
    u_hat: (3, n, n, n) complex
    """
    n = u_hat.shape[-1]
    kx, ky, kz = wave_number_grids(n)
    k2 = kx * kx + ky * ky + kz * kz
    safe = k2.copy()
    safe[0, 0, 0] = 1.0
    div = kx * u_hat[0] + ky * u_hat[1] + kz * u_hat[2]
    factor = div / safe
    out = np.empty_like(u_hat)
    out[0] = u_hat[0] - kx * factor
    out[1] = u_hat[1] - ky * factor
    out[2] = u_hat[2] - kz * factor
    out[:, 0, 0, 0] = u_hat[:, 0, 0, 0]
    return out


def divergence_hat(u_hat: np.ndarray) -> np.ndarray:
    n = u_hat.shape[-1]
    kx, ky, kz = wave_number_grids(n)
    return kx * u_hat[0] + ky * u_hat[1] + kz * u_hat[2]


def divergence_l2(u_hat: np.ndarray) -> float:
    d = divergence_hat(u_hat)
    return float(np.sqrt(np.sum(np.abs(d) ** 2)))
