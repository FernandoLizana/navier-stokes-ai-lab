"""
Coefficient-space evolution strategy (CMA-lite) to maximize enstrophy(T)
on the fixed-energy divergence-free sphere — attack on C-0001.

Pure NumPy (μ/λ)-ES with rank-μ covariance adaptation (simplified).
Evidence ≤ N2. Not a maximizer certificate. Not Clay.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ns_exploration.optimization.adjoint_ascent import enstrophy_integrator
from ns_exploration.spectral.fourier_conventions import enforce_reality, kinetic_energy_from_hat
from ns_exploration.spectral.leray import leray_project_hat


def _pack_divfree_low_modes(n: int, k_max: int = 3) -> list[tuple[int, int, int, int]]:
    """
    Degrees of freedom: for each k with 0 < |k|_∞ ≤ k_max and k in FFT ordering,
    store 2 real params per free tangential direction (2 per mode for div-free in 3D).
    We use a simpler encoding: full complex vector then Leray+reality+energy project.
    Returns list of (component, ix, iy, iz) for independent complex coeffs with
    lexicographic k and not the conjugate partner.
    """
    from ns_exploration.spectral.fourier_conventions import wave_numbers_1d

    kx = wave_numbers_1d(n)
    dofs: list[tuple[int, int, int, int]] = []
    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                k = (kx[ix], kx[iy], kx[iz])
                if k == (0.0, 0.0, 0.0):
                    continue
                if max(abs(k[0]), abs(k[1]), abs(k[2])) > k_max + 1e-9:
                    continue
                # Keep half-space to ease reality: first nonzero of (ix,iy,iz) vs conjugate
                # Use: flatten index ordering — keep if (ix,iy,iz) <= conjugate indices lex
                jx = (n - ix) % n
                jy = (n - iy) % n
                jz = (n - iz) % n
                if (ix, iy, iz) > (jx, jy, jz):
                    continue
                if (ix, iy, iz) == (jx, jy, jz):
                    # Nyquist / self-conjugate: real only — still allow 3 comps then Leray
                    for c in range(3):
                        dofs.append((c, ix, iy, iz))
                else:
                    for c in range(3):
                        dofs.append((c, ix, iy, iz))
    return dofs


def theta_to_uhat(
    theta: np.ndarray,
    n: int,
    dofs: list[tuple[int, int, int, int]],
    energy_target: float,
) -> np.ndarray:
    """Map real parameter vector → div-free real field with fixed energy."""
    u = np.zeros((3, n, n, n), dtype=np.complex128)
    # theta layout: for each dof, (re, im) except we always use 2 floats per dof entry
    need = 2 * len(dofs)
    if theta.size < need:
        raise ValueError(f"theta length {theta.size} < {need}")
    for i, (c, ix, iy, iz) in enumerate(dofs):
        re = theta[2 * i]
        im = theta[2 * i + 1]
        u[c, ix, iy, iz] += re + 1j * im
    u = enforce_reality(u)
    u = leray_project_hat(u)
    e = kinetic_energy_from_hat(u)
    if e <= 1e-30:
        # fallback tiny random
        u = np.zeros_like(u)
        u[0, 1, 0, 0] = 0.5
        u[0, -1, 0, 0] = 0.5
        u = leray_project_hat(u)
        e = kinetic_energy_from_hat(u)
    return u * np.sqrt(energy_target / e)


@dataclass
class CMAESResult:
    n: int
    k_max: int
    generations: int
    best_J: float
    best_generation: int
    history_best: list[float]
    n_evals: int
    exceeded_M: bool
    M_target: float
    evidence_level: str = "N2"
    route: str = "B"
    notes: str = "Simplified (μ/λ)-ES; not full Hansen CMA-ES."

    def as_dict(self) -> dict:
        return asdict(self)


def evolve_enstrophy_es(
    n: int = 12,
    k_max: int = 3,
    energy_target: float = 0.5,
    nu: float = 0.1,
    dt: float = 1e-3,
    t_end: float = 0.02,
    generations: int = 25,
    population: int = 12,
    mu: int = 4,
    sigma0: float = 0.15,
    seed: int = 0,
    M_target: float = 10.994055804604205,
    integrator: str = "etd_rk2",
) -> tuple[np.ndarray, CMAESResult]:
    """
    Maximize enstrophy(T) over low-mode coefficient sphere.
    """
    rng = np.random.default_rng(seed)
    dofs = _pack_divfree_low_modes(n, k_max=k_max)
    dim = 2 * len(dofs)
    mean = rng.normal(scale=0.05, size=dim)
    sigma = sigma0
    # Diagonal covariance adaptation (cheap CMA-lite)
    diag_c = np.ones(dim)

    def fitness(theta: np.ndarray) -> float:
        uh = theta_to_uhat(theta, n, dofs, energy_target)
        return enstrophy_integrator(uh, nu, dt, t_end, integrator)

    best_J = -np.inf
    best_theta = mean.copy()
    best_gen = 0
    history: list[float] = []
    n_evals = 0

    for g in range(generations):
        thetas = []
        scores = []
        for _ in range(population):
            z = rng.normal(size=dim) * np.sqrt(diag_c)
            theta = mean + sigma * z
            j = fitness(theta)
            n_evals += 1
            thetas.append(theta)
            scores.append(j)
            if j > best_J:
                best_J = j
                best_theta = theta.copy()
                best_gen = g
        history.append(best_J)
        # Select top mu
        order = np.argsort(scores)[::-1]
        selected = [thetas[i] for i in order[:mu]]
        weights = np.array([math_log_mu(mu, i) for i in range(mu)])
        weights = weights / weights.sum()
        new_mean = sum(w * th for w, th in zip(weights, selected))
        # Update diagonal covariance toward selected steps
        for w, th in zip(weights, selected):
            step = (th - mean) / (sigma + 1e-30)
            diag_c = (1 - w * 0.2) * diag_c + w * 0.2 * (step**2)
        diag_c = np.clip(diag_c, 1e-4, 1e4)
        mean = new_mean
        # Adapt sigma: expand if best improving, else shrink mildly
        if g > 0 and history[-1] > history[-2] + 1e-6:
            sigma *= 1.05
        else:
            sigma *= 0.96
        sigma = float(np.clip(sigma, 1e-3, 1.0))
        if best_J > M_target:
            break

    best_u = theta_to_uhat(best_theta, n, dofs, energy_target)
    result = CMAESResult(
        n=n,
        k_max=k_max,
        generations=len(history),
        best_J=float(best_J),
        best_generation=best_gen,
        history_best=history,
        n_evals=n_evals,
        exceeded_M=bool(best_J > M_target),
        M_target=M_target,
    )
    return best_u, result


def math_log_mu(mu: int, i: int) -> float:
    """Positive decreasing weights ~ log(μ+1)-log(i+1)."""
    return float(np.log(mu + 0.5) - np.log(i + 1))
