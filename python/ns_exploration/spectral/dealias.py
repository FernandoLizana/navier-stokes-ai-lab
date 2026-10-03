"""Orszag 2/3 dealiasing rule."""

from __future__ import annotations

import numpy as np

from ns_exploration.spectral.fourier_conventions import wave_number_grids


def dealias_mask(n: int) -> np.ndarray:
    """True where mode is kept (|k_i| < n/3 for all i)."""
    kx, ky, kz = wave_number_grids(n)
    kmax = n / 3.0
    return (np.abs(kx) < kmax) & (np.abs(ky) < kmax) & (np.abs(kz) < kmax)


def apply_dealias(u_hat: np.ndarray, mask: np.ndarray | None = None) -> np.ndarray:
    if mask is None:
        mask = dealias_mask(u_hat.shape[-1])
    out = u_hat.copy()
    out[:, ~mask] = 0.0
    return out
