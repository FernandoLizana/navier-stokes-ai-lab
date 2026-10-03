"""
NS-MRL Fourier conventions — SINGLE SOURCE OF TRUTH.

Domain: T^3 = [0, 2π)^3 with integer wavevectors k ∈ ℤ³.
Representation: u(x) = Σ_k û_k exp(i k · x)

Discrete N³ grid, N even recommended.
Coefficient array û = fftn(u_phys) / N³  (true Fourier coefficients).
Reconstruction: u_phys = ifftn(û * N³).

Parseval: (1/(2π)³) ∫ |u|² dx = Σ_k |û_k|²
Kinetic energy E = (1/2) Σ_k |û_k|² = (1/2) mean(|u_phys|²)

All spectral modules MUST import from here. Dual normalizations are forbidden.
Route: B (periodic). Evidence of numerics using this: ≤ N2 unless validated.
"""

from __future__ import annotations

import numpy as np

TWOPI = 2.0 * np.pi
DOMAIN_PERIOD = TWOPI
NORMALIZATION = "true_fourier_coefficients_T3_0_2pi"


def wave_numbers_1d(n: int) -> np.ndarray:
    """Integer wave numbers for an n-point periodic grid on [0, 2π)."""
    return np.fft.fftfreq(n, d=1.0 / n).astype(np.float64)


def wave_number_grids(n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    kx = wave_numbers_1d(n)
    ky = wave_numbers_1d(n)
    kz = wave_numbers_1d(n)
    return np.meshgrid(kx, ky, kz, indexing="ij")


def k_squared(n: int) -> np.ndarray:
    kx, ky, kz = wave_number_grids(n)
    return kx * kx + ky * ky + kz * kz


def fft(u_phys: np.ndarray) -> np.ndarray:
    """Physical → Fourier coefficients. u_phys shape (n, n, n)."""
    n = u_phys.shape[-1]
    return np.fft.fftn(u_phys, axes=(-3, -2, -1)) / (n**3)


def ifft(u_hat: np.ndarray) -> np.ndarray:
    """Fourier coefficients → physical. Returns real part for real fields."""
    n = u_hat.shape[-1]
    return np.fft.ifftn(u_hat * (n**3), axes=(-3, -2, -1)).real


def fft_vector(u_phys: np.ndarray) -> np.ndarray:
    """u_phys shape (3, n, n, n) → û shape (3, n, n, n)."""
    return np.stack([fft(u_phys[i]) for i in range(3)], axis=0)


def ifft_vector(u_hat: np.ndarray) -> np.ndarray:
    return np.stack([ifft(u_hat[i]) for i in range(3)], axis=0)


def enforce_reality(u_hat: np.ndarray) -> np.ndarray:
    """Project to conjugate-symmetric coefficients of a real field via round-trip."""
    return fft_vector(ifft_vector(u_hat))


def kinetic_energy_from_hat(u_hat: np.ndarray) -> float:
    """E = (1/2) Σ |û_k|² (sum over components and modes)."""
    return 0.5 * float(np.sum(np.abs(u_hat) ** 2))


def kinetic_energy_from_phys(u_phys: np.ndarray) -> float:
    return 0.5 * float(np.mean(np.sum(u_phys**2, axis=0)))
