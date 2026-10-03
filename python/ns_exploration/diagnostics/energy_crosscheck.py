"""
Independent energy-balance cross-checks (exploratory ≤ N3).

1) Physical-space mean energy vs spectral Parseval.
2) 2D NS embedded in 3D (w≡0, ∂z≡0): enstrophy/energy evolution
   compared between full 3D spectral operator and a reduced 2D spectral kernel.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np

from ns_exploration.spectral.fourier_conventions import (
    fft,
    fft_vector,
    ifft,
    ifft_vector,
    kinetic_energy_from_hat,
    kinetic_energy_from_phys,
    wave_number_grids,
)
from ns_exploration.spectral.integrators import step_rk4
from ns_exploration.spectral.leray import divergence_l2, leray_project_hat
from ns_exploration.spectral.operators import curl_hat, nonlinear_hat


@dataclass
class EnergyCrossReport:
    name: str
    n: int
    energy_spectral: float
    energy_physical: float
    rel_diff: float
    div_l2: float
    notes: str
    evidence_level: str = "N2"
    route: str = "B"

    def as_dict(self) -> dict:
        return asdict(self)


def parseval_crosscheck(u_hat: np.ndarray) -> EnergyCrossReport:
    e_s = kinetic_energy_from_hat(u_hat)
    e_p = kinetic_energy_from_phys(ifft_vector(u_hat))
    rel = abs(e_s - e_p) / (abs(e_s) + 1e-30)
    return EnergyCrossReport(
        name="parseval_spectral_vs_physical",
        n=u_hat.shape[-1],
        energy_spectral=e_s,
        energy_physical=e_p,
        rel_diff=rel,
        div_l2=divergence_l2(u_hat),
        notes="Must agree up to roundoff under NS-MRL Fourier convention.",
    )


def embedded_2d_taylor_green(n: int, a: float = 1.0) -> np.ndarray:
    """2D TG vortex extruded in z (independent of z), w=0, div-free."""
    from ns_exploration.spectral.fourier_conventions import TWOPI

    x = np.linspace(0.0, TWOPI, n, endpoint=False)
    X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
    u = a * np.sin(X) * np.cos(Y)
    v = -a * np.cos(X) * np.sin(Y)
    w = np.zeros_like(u)
    return leray_project_hat(fft_vector(np.stack([u, v, w], axis=0)))


def nonlinear_2d_reduced(u_hat: np.ndarray, dealias: bool = True) -> np.ndarray:
    """
    Reduced 2D convective operator on (u,v) ignoring z-derivatives,
    then embedded back as (Nu, Nv, 0) with Leray in 3D.
    Independent code path from operators.nonlinear_hat (no z derivatives used).
    """
    from ns_exploration.spectral.dealias import apply_dealias, dealias_mask

    n = u_hat.shape[-1]
    kx, ky, _kz = wave_number_grids(n)
    u = ifft(u_hat[0])
    v = ifft(u_hat[1])
    ux = ifft(1j * kx * u_hat[0])
    uy = ifft(1j * ky * u_hat[0])
    vx = ifft(1j * kx * u_hat[1])
    vy = ifft(1j * ky * u_hat[1])
    conv_u = u * ux + v * uy
    conv_v = u * vx + v * vy
    conv_hat = np.stack(
        [
            fft(conv_u),
            fft(conv_v),
            np.zeros((n, n, n), dtype=np.complex128),
        ],
        axis=0,
    )
    if dealias:
        conv_hat = apply_dealias(conv_hat, dealias_mask(n))
    return -leray_project_hat(conv_hat)


def embedded_2d_operator_agreement(n: int = 24) -> dict:
    """Compare full 3D nonlinear vs reduced 2D path on z-invariant field."""
    uh = embedded_2d_taylor_green(n)
    N3 = nonlinear_hat(uh, dealias=True)
    N2 = nonlinear_2d_reduced(uh, dealias=True)
    diff = float(np.sqrt(np.sum(np.abs(N3 - N2) ** 2)))
    scale = float(np.sqrt(np.sum(np.abs(N3) ** 2))) + 1e-30
    # Energy channel: Re <u, N> should be ~0
    inner3 = float(np.sum((uh.conj() * N3).real))
    inner2 = float(np.sum((uh.conj() * N2).real))
    return {
        "n": n,
        "operator_rel_diff": diff / scale,
        "energy_channel_3d": inner3,
        "energy_channel_2d": inner2,
        "div_l2": divergence_l2(uh),
        "evidence_level": "N2",
        "route": "B",
        "notes": "2D embedded subclass only — does not imply 3D regularity.",
    }


def energy_decay_crosscheck(n: int = 16, nu: float = 0.05, dt: float = 1e-3, steps: int = 40) -> dict:
    """
    Compare ΔE / Δt to -ν ||∇u||² ≈ -2ν enstrophy using spectral diagnostics.
    Independent of integrator internals beyond using RK4 steps.
    """
    uh = embedded_2d_taylor_green(n)
    e0 = kinetic_energy_from_hat(uh)
    w0 = curl_hat(uh)
    ens0 = 0.5 * float(np.sum(np.abs(w0) ** 2))
    for _ in range(steps):
        uh = step_rk4(uh, dt, nu, dealias=True)
    e1 = kinetic_energy_from_hat(uh)
    # Trapezoid estimate of dissipation from endpoints is crude; use mean of start ens
    # Better: accumulate instantaneous balance each step
    uh = embedded_2d_taylor_green(n)
    dissip_trap = 0.0
    e_prev = kinetic_energy_from_hat(uh)
    energy_residuals = []
    for _ in range(steps):
        w = curl_hat(uh)
        ens = 0.5 * float(np.sum(np.abs(w) ** 2))
        pred = e_prev - dt * (2.0 * nu * ens)  # forward Euler predictor of identity
        uh = step_rk4(uh, dt, nu, dealias=True)
        e_new = kinetic_energy_from_hat(uh)
        energy_residuals.append(e_new - pred)
        dissip_trap += 2.0 * nu * ens * dt
        e_prev = e_new
    e1 = kinetic_energy_from_hat(uh)
    return {
        "n": n,
        "delta_E": e1 - e0,
        "minus_integrated_dissipation_proxy": -dissip_trap,
        "mean_abs_step_residual": float(np.mean(np.abs(energy_residuals))),
        "max_abs_step_residual": float(np.max(np.abs(energy_residuals))),
        "ens0": ens0,
        "evidence_level": "N2",
        "route": "B",
        "notes": "Step residual mixes time-discretization error with dealiasing; not a proof.",
    }
