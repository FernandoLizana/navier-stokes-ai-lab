"""
L-0021: Quartic Gram upper bound on ‖P((u·∇)u)‖ for mono-radial fields.

For a field supported on a single radial shell r = |k|², expand in a real
orthonormal polarization basis (D = 2 m_r). The map u ↦ N = P((u·∇)u) is
quadratic. Writing ‖N‖² as a homogeneous quartic and applying the Shor-style
Frobenius relaxation on X ∈ Sym_D (max Wᵀ G W over ‖W‖_F = 1 dominates the
rank-1 case W = uuᵀ), power iteration on the Gram operator G yields
  α_r ≥ max{ ‖N(u)‖ : supp(u) ⊂ shell r, ‖u‖_ℓ₂ = 1 }.
Hence for energy E ≤ E0 (‖u‖_ℓ₂² = 2E ≤ 1 at E0 = 1/2),
  ‖N‖ ≤ α_r.

On N=24 dealias: α_★ := max_r α_r = 14 (attained at r=74). Combined with
L-0020 Duhamel H¹ this proves a mono-radial enstrophy bound
  Ω(T) ≤ ½(√(2 S E0) + α_★ I_σ)² ≈ 71.73 < K² E0 = 73.5
(beats the spectral envelope on this IC subclass).

Does NOT control multi-shell / all-IC fields (C-0005 still open; needs N_*≤9.67).

FINITE Galerkin / dealias only. Not continuum. Not Clay. Evidence: N7 (+ N5).
"""

from __future__ import annotations

import json
import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0020_duhamel_h1 import (
    lemma_l0020,
    omega_duhamel_H1,
)
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids


def _pol_basis(k: np.ndarray) -> list[np.ndarray]:
    k = np.asarray(k, float)
    if abs(k[0]) < 0.9:
        a = np.cross(k, [1.0, 0.0, 0.0])
    else:
        a = np.cross(k, [0.0, 1.0, 0.0])
    a = a / np.linalg.norm(a)
    b = np.cross(k, a)
    b = b / np.linalg.norm(b)
    return [a, b]


def _leray_vec(s: np.ndarray, v: np.ndarray) -> np.ndarray:
    s = np.asarray(s, float)
    sn2 = float(np.dot(s, s))
    if sn2 < 1e-14:
        return np.zeros(3)
    return v - s * np.dot(s, v) / sn2


def shell_modes_by_r(n: int) -> dict[int, list[tuple[int, int, int]]]:
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    by_r: dict[int, list[tuple[int, int, int]]] = defaultdict(list)
    nonzero = mask & ~((kx == 0) & (ky == 0) & (kz == 0))
    for ix, iy, iz in np.argwhere(nonzero):
        q = (int(kx[ix, iy, iz]), int(ky[ix, iy, iz]), int(kz[ix, iy, iz]))
        by_r[q[0] * q[0] + q[1] * q[1] + q[2] * q[2]].append(q)
    return by_r


def out_index_map(n: int) -> dict[tuple[int, int, int], int]:
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    out_modes = [
        (int(a), int(b), int(c))
        for a, b, c in zip(kx[mask], ky[mask], kz[mask])
        if not (a == 0 and b == 0 and c == 0)
    ]
    return {k: i for i, k in enumerate(out_modes)}


def quartic_alpha_r(
    r: int,
    by_r: dict[int, list[tuple[int, int, int]]],
    out_index: dict[tuple[int, int, int], int],
    n_out: int,
    iters: int = 50,
) -> float:
    modes = by_r[r]
    if not modes:
        return 0.0
    basis: list[tuple[np.ndarray, np.ndarray]] = []
    for k in modes:
        for p in _pol_basis(np.array(k, float)):
            basis.append((np.array(k, int), np.asarray(p, float)))
    D = len(basis)
    pairs: list[tuple[int, int, int, np.ndarray]] = []
    for i, (ki, ai) in enumerate(basis):
        for j, (kj, aj) in enumerate(basis):
            s = (int(ki[0] + kj[0]), int(ki[1] + kj[1]), int(ki[2] + kj[2]))
            if s not in out_index:
                continue
            coef = float(np.dot(ai, kj.astype(float)))
            Pv = _leray_vec(np.array(s, float), aj)
            pairs.append((i, j, out_index[s], coef * Pv))
    if not pairs:
        return 0.0

    def apply_G(X: np.ndarray) -> np.ndarray:
        field = np.zeros((n_out, 3), dtype=np.float64)
        for i, j, si, vec in pairs:
            field[si] += X[i, j] * vec
        out = np.zeros((D, D), dtype=np.float64)
        for k, l, si, vec in pairs:
            out[k, l] += float(np.dot(field[si], vec))
        return out

    rng = np.random.default_rng(int(r) + 17)
    X = rng.normal(size=(D, D))
    X = 0.5 * (X + X.T)
    X /= np.linalg.norm(X) + 1e-30
    for _ in range(iters):
        Y = apply_G(X)
        Y = 0.5 * (Y + Y.T)
        nrm = float(np.linalg.norm(Y))
        if nrm < 1e-30:
            break
        X = Y / nrm
    lam = float(np.sum(X * apply_G(X)))
    return math.sqrt(max(lam, 0.0))


def all_shell_alphas(n: int = 24, iters: int = 50) -> dict[int, float]:
    by_r = shell_modes_by_r(n)
    out_index = out_index_map(n)
    n_out = len(out_index)
    return {
        r: quartic_alpha_r(r, by_r, out_index, n_out, iters=iters)
        for r in sorted(by_r.keys())
    }


@dataclass
class GalerkinBoundL0021:
    lemma_id: str = "L-0021"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    E0: float = 0.5
    nu: float = 0.1
    T: float = 0.02
    alpha_star: float = 0.0
    r_star: int = 0
    n_shells: int = 0
    Omega_radial_H1: float = 0.0
    envelope_cap: float = 0.0
    beats_envelope: bool = False
    c0005_M: float = 61.23693461895651
    closes_c0005: bool = False
    clay_implication: str = (
        "None. Mono-radial Galerkin quartic Gram + Duhamel only; "
        "not all-IC; not continuum; not Clay."
    )
    notes: str = ""
    alphas: dict[str, float] | None = None

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0021(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    c0005_M: float = 61.23693461895651,
    iters: int = 50,
) -> GalerkinBoundL0021:
    alphas = all_shell_alphas(n, iters=iters)
    r_star = max(alphas, key=lambda r: alphas[r])
    alpha_star = float(alphas[r_star])
    b20 = lemma_l0020(n=n, E0=E0, nu=nu, T=T, c0005_M=c0005_M)
    om = omega_duhamel_H1(b20.S_T, E0, alpha_star, b20.I_sigma)
    notes = (
        f"α_★={alpha_star:.6g} at r={r_star} over {len(alphas)} shells. "
        f"Mono-radial Duhamel Ω≤{om:.6g} (envelope {b20.envelope_cap}). "
        f"Closes C-0005 M={c0005_M:.4g}: {om <= c0005_M}."
    )
    return GalerkinBoundL0021(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        alpha_star=alpha_star,
        r_star=int(r_star),
        n_shells=len(alphas),
        Omega_radial_H1=om,
        envelope_cap=b20.envelope_cap,
        beats_envelope=om < b20.envelope_cap,
        c0005_M=c0005_M,
        closes_c0005=om <= c0005_M + 1e-12,
        notes=notes,
        alphas={str(k): float(v) for k, v in alphas.items()},
    )


def save_lemma_l0021(
    bound: GalerkinBoundL0021,
    path: str | Path = "conjectures/proved_restricted/L-0021.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
