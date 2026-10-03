"""
Build a single |k|²=K²_max dealias mode (Stokes IC) used to refute C-0002.

For a lone Fourier mode, P((u·∇)u)=0, so the Galerkin NS flow is exact Stokes:
  Ω(t) = K² E0 exp(-2 ν K² t).
At N=16, K²=75, E0=0.5, ν=0.1, T=0.02: Ω(T)≈27.7807 > M_C0002≈20.32.
"""

from __future__ import annotations

import math

import numpy as np

from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import (
    enforce_reality,
    k_squared,
    kinetic_energy_from_hat,
    wave_number_grids,
)
from ns_exploration.spectral.leray import leray_project_hat
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def stokes_max_mode_ic(
    n: int = 16,
    E0: float = 0.5,
) -> tuple[np.ndarray, dict]:
    """Unit-energy-scaled real Leray field supported on one ±k pair with |k|² = K²_max."""
    K2, _ = full_dealias_exact_stats(n)
    mask = dealias_mask(n)
    k2 = k_squared(n)
    kx, ky, kz = wave_number_grids(n)
    idxs = np.argwhere((np.rint(k2) == K2) & mask)
    if idxs.size == 0:
        raise RuntimeError(f"no dealias mode with |k|²={K2} on N={n}")
    ix, iy, iz = map(int, idxs[0])
    kvec = np.array([kx[ix, iy, iz], ky[ix, iy, iz], kz[ix, iy, iz]], dtype=float)
    if abs(kvec[0]) < 0.9:
        a = np.cross(kvec, np.array([1.0, 0.0, 0.0]))
    else:
        a = np.cross(kvec, np.array([0.0, 1.0, 0.0]))
    a /= np.linalg.norm(a)

    u = np.zeros((3, n, n, n), dtype=np.complex128)
    c = 1.0 / math.sqrt(2.0)
    u[:, ix, iy, iz] = c * a
    neg = np.argwhere(
        np.isclose(kx, -kvec[0]) & np.isclose(ky, -kvec[1]) & np.isclose(kz, -kvec[2])
    )[0]
    u[:, int(neg[0]), int(neg[1]), int(neg[2])] = c * a
    u = enforce_reality(leray_project_hat(u))
    e = kinetic_energy_from_hat(u)
    u *= math.sqrt(E0 / e)

    meta = {
        "n": n,
        "K2": int(K2),
        "k": [float(kvec[0]), float(kvec[1]), float(kvec[2])],
        "E0": float(kinetic_energy_from_hat(u)),
        "Omega0": float(K2) * E0,
        "formula_Omega_T": "K2*E0*exp(-2*nu*K2*T)",
    }
    return u, meta


def stokes_Omega_exact(K2: int, E0: float, nu: float, t: float) -> float:
    return float(K2) * E0 * math.exp(-2.0 * nu * float(K2) * t)
