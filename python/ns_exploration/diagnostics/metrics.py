"""Exploratory diagnostics. Not singularity proofs."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ns_exploration.spectral.fourier_conventions import (
    ifft_vector,
    kinetic_energy_from_hat,
    wave_number_grids,
)
from ns_exploration.spectral.leray import divergence_l2
from ns_exploration.spectral.operators import curl_hat


@dataclass
class Diagnostics:
    energy: float
    enstrophy: float
    palinstrophy: float
    helicity: float
    omega_inf: float
    u_inf: float
    dissipation: float
    div_l2: float
    spectral_tail_energy: float
    analyticity_radius_proxy: float

    def as_dict(self) -> dict:
        return asdict(self)


def compute_diagnostics(u_hat: np.ndarray, nu: float) -> Diagnostics:
    n = u_hat.shape[-1]
    energy = kinetic_energy_from_hat(u_hat)
    w_hat = curl_hat(u_hat)
    enstrophy = 0.5 * float(np.sum(np.abs(w_hat) ** 2))
    # palinstrophy ~ 0.5 ||∇ω||²
    kx, ky, kz = wave_number_grids(n)
    k2 = kx * kx + ky * ky + kz * kz
    pal = 0.5 * float(np.sum(k2 * np.sum(np.abs(w_hat) ** 2, axis=0)))
    u = ifft_vector(u_hat)
    w = ifft_vector(w_hat)
    helicity = float(np.mean(np.sum(u * w, axis=0)))
    omega_inf = float(np.max(np.sqrt(np.sum(w**2, axis=0))))
    u_inf = float(np.max(np.sqrt(np.sum(u**2, axis=0))))
    dissipation = 2.0 * nu * enstrophy  # since dE/dt = -ν ||∇u||² = -2ν enstrophy for our E
    # spectral tail: energy in outer third of modes
    kmax = n // 2
    kr = np.sqrt(k2)
    tail = kr > (2.0 * kmax / 3.0)
    tail_e = 0.5 * float(np.sum(np.abs(u_hat[:, tail]) ** 2))
    # crude analyticity proxy: fit log|E(k)| vs k on shells — use highest shell with energy
    shells = {}
    kint = np.rint(kr).astype(int)
    for kk in range(1, kmax):
        m = kint == kk
        if np.any(m):
            shells[kk] = 0.5 * float(np.sum(np.abs(u_hat[:, m]) ** 2))
    rho = 0.0
    if len(shells) >= 3:
        ks = np.array(sorted(shells.keys()), dtype=float)
        es = np.array([shells[int(k)] for k in ks]) + 1e-30
        # linear fit log E = a - ρ k
        coeff = np.polyfit(ks[es > 0], np.log(es[es > 0]), 1)
        rho = float(max(0.0, -coeff[0]))
    return Diagnostics(
        energy=energy,
        enstrophy=enstrophy,
        palinstrophy=pal,
        helicity=helicity,
        omega_inf=omega_inf,
        u_inf=u_inf,
        dissipation=dissipation,
        div_l2=divergence_l2(u_hat),
        spectral_tail_energy=tail_e,
        analyticity_radius_proxy=rho,
    )
