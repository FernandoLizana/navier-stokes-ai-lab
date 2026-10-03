"""Reproducible divergence-free smooth initial conditions on T³. Route B."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass

import numpy as np

from ns_exploration.spectral.fourier_conventions import (
    TWOPI,
    fft_vector,
    kinetic_energy_from_hat,
)
from ns_exploration.spectral.leray import leray_project_hat
from ns_exploration.spectral.operators import curl_hat


@dataclass
class ICMetadata:
    name: str
    n: int
    seed: int | None
    energy: float
    enstrophy: float
    helicity: float
    symmetries: list[str]
    resolution: int
    content_hash: str
    license: str = "MIT"
    provenance: str = "ns_exploration.initial_conditions"
    route: str = "B"


def _grid(n: int):
    x = np.linspace(0.0, TWOPI, n, endpoint=False)
    return np.meshgrid(x, x, x, indexing="ij")


def _meta(name: str, u_hat: np.ndarray, seed: int | None, symmetries: list[str]) -> ICMetadata:
    w_hat = curl_hat(u_hat)
    payload = np.concatenate([u_hat.real.ravel(), u_hat.imag.ravel()])
    h = hashlib.sha256(payload.tobytes()).hexdigest()[:16]
    return ICMetadata(
        name=name,
        n=u_hat.shape[-1],
        seed=seed,
        energy=kinetic_energy_from_hat(u_hat),
        enstrophy=0.5 * float(np.sum(np.abs(w_hat) ** 2)),
        helicity=0.0,  # filled by caller if needed
        symmetries=symmetries,
        resolution=u_hat.shape[-1],
        content_hash=h,
    )


def taylor_green(n: int, a: float = 1.0) -> tuple[np.ndarray, ICMetadata]:
    """Classic Taylor–Green vortex (divergence-free)."""
    X, Y, Z = _grid(n)
    u = a * np.sin(X) * np.cos(Y) * np.cos(Z)
    v = -a * np.cos(X) * np.sin(Y) * np.cos(Z)
    w = np.zeros_like(u)
    u_hat = leray_project_hat(fft_vector(np.stack([u, v, w], axis=0)))
    meta = _meta("taylor_green", u_hat, None, ["TG"])
    return u_hat, meta


def abc_flow(n: int, A: float = 1.0, B: float = 1.0, C: float = 1.0) -> tuple[np.ndarray, ICMetadata]:
    """Arnold–Beltrami–Childress Beltrami field (steady Euler when ν=0,f=0)."""
    X, Y, Z = _grid(n)
    u = A * np.sin(Z) + C * np.cos(Y)
    v = B * np.sin(X) + A * np.cos(Z)
    w = C * np.sin(Y) + B * np.cos(X)
    u_hat = leray_project_hat(fft_vector(np.stack([u, v, w], axis=0)))
    meta = _meta("abc", u_hat, None, ["ABC", "Beltrami"])
    return u_hat, meta


def random_div_free(
    n: int,
    seed: int = 0,
    k_peak: int = 3,
    energy_target: float = 0.5,
) -> tuple[np.ndarray, ICMetadata]:
    """Gaussian random Fourier field, projected divergence-free, scaled to energy."""
    rng = np.random.default_rng(seed)
    shape = (3, n, n, n)
    real = rng.normal(size=shape)
    imag = rng.normal(size=shape)
    u_hat = real + 1j * imag
    # Spectral envelope
    from ns_exploration.spectral.fourier_conventions import k_squared

    k2 = k_squared(n)
    env = np.exp(-0.5 * (np.sqrt(k2) - k_peak) ** 2)
    u_hat *= env
    u_hat[:, 0, 0, 0] = 0.0
    u_hat = leray_project_hat(u_hat)
    # Reality via round-trip
    from ns_exploration.spectral.fourier_conventions import enforce_reality

    u_hat = enforce_reality(u_hat)
    u_hat = leray_project_hat(u_hat)
    e = kinetic_energy_from_hat(u_hat)
    if e > 0:
        u_hat *= np.sqrt(energy_target / e)
    meta = _meta("random_div_free", u_hat, seed, ["random"])
    meta.energy = kinetic_energy_from_hat(u_hat)
    return u_hat, meta


def vortex_tubes_periodic(
    n: int,
    amplitude: float = 1.0,
    sigma: float = 0.5,
    n_tubes: int = 2,
) -> tuple[np.ndarray, ICMetadata]:
    """
    Smooth periodic anti-parallel vortex tubes (Gaussian cores along z), then Leray project.
    Result is smooth and divergence-free (projection may slightly alter tubes).
    """
    X, Y, Z = _grid(n)
    # Streamfunction-like construction: ω ≈ (0,0,ωz), u from Biot–Savart spectral
    omega_z = np.zeros((n, n, n))
    centers = [(np.pi * 0.5, np.pi), (np.pi * 1.5, np.pi)][:n_tubes]
    signs = [1.0, -1.0]
    for (cx, cy), s in zip(centers, signs):
        # periodic minimum-image distance in x,y
        dx = np.angle(np.exp(1j * (X - cx)))
        dy = np.angle(np.exp(1j * (Y - cy)))
        omega_z += s * amplitude * np.exp(-(dx**2 + dy**2) / (2 * sigma**2))
    # Build vorticity hat and recover velocity: û = ik×ω̂ / |k|² (divergence-free)
    from ns_exploration.spectral.fourier_conventions import fft, wave_number_grids

    wz_hat = fft(omega_z)
    kx, ky, kz = wave_number_grids(n)
    k2 = kx * kx + ky * ky + kz * kz
    safe = k2.copy()
    safe[0, 0, 0] = 1.0
    # ω = (0,0,ωz) → u = curl(ψ) with -Δψ=ω ⇒ û = ik × ω̂ / k²
    ox, oy, oz = 0j, 0j, wz_hat
    ux = 1j * (ky * oz - kz * oy) / safe
    uy = 1j * (kz * ox - kx * oz) / safe
    uz = 1j * (kx * oy - ky * ox) / safe
    u_hat = np.stack([ux, uy, uz], axis=0)
    u_hat[:, 0, 0, 0] = 0.0
    u_hat = leray_project_hat(u_hat)
    from ns_exploration.spectral.fourier_conventions import enforce_reality

    # Reality round-trip can reintroduce O(eps) divergence — re-project.
    u_hat = leray_project_hat(enforce_reality(u_hat))
    meta = _meta("vortex_tubes_periodic", u_hat, None, ["tubes"])
    return u_hat, meta


def save_ic(path: str, u_hat: np.ndarray, meta: ICMetadata) -> None:
    from pathlib import Path

    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        u_hat_real=u_hat.real,
        u_hat_imag=u_hat.imag,
        meta=json.dumps(asdict(meta)),
    )
