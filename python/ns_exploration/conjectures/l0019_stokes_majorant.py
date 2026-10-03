"""
L-0019: Stokes spectral majorant (linear) + floor for C-0005.

On the dealias Galerkin space, the Stokes semigroup gives
  û_k(t) = e^{-ν |k|² t} û_k(0)
for each mode when the nonlinear term vanishes (e.g. single Fourier mode).
Hence
  Ω(t) = (1/2) Σ |k|² |û_k(t)|² ≤ S(t) E(0),
  S(t) := max_{k ∈ mask} |k|² exp(-2 ν |k|² t).

For N=24, ν=0.1, T=0.02: the max is attained at K²=147 and
  S(T) E0 ≈ 40.8246 =: Stokes_floor,
which is the greatest lower bound any all-IC claim at that T must clear
(and is realized by the max-mode IC).

This lemma proves the *linear* majorant only. It does NOT bound full NS
(nonlinear stretching may exceed the Stokes curve). Used to:
  - refute all-IC conjectures with M < Stokes_floor,
  - place C-0005 (M≈61.24) strictly above the floor and below the L-0018 envelope.

FINITE Galerkin / dealias only. Not continuum. Not Clay. Evidence: N7 (Stokes).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def stokes_S(n: int, nu: float, t: float) -> tuple[float, int]:
    """Return (S(t), argmax |k|²) on the dealias mask."""
    mask = dealias_mask(n)
    k2 = k_squared(n)
    vals = np.rint(k2[mask]).astype(np.int64)
    if vals.size == 0:
        return 0.0, 0
    best_S = -1.0
    best_r = 0
    for r in np.unique(vals):
        r = int(r)
        if r <= 0:
            continue
        S = float(r) * math.exp(-2.0 * nu * float(r) * t)
        if S > best_S:
            best_S = S
            best_r = r
    # Cross-check against K²_max (for t small / νT small, max is at K²_max)
    K2, _ = full_dealias_exact_stats(n)
    S_K = float(K2) * math.exp(-2.0 * nu * float(K2) * t)
    if S_K >= best_S - 1e-15:
        return S_K, int(K2)
    return best_S, best_r


def stokes_floor(n: int, E0: float, nu: float, t: float) -> float:
    S, _ = stokes_S(n, nu, t)
    return S * E0


@dataclass
class GalerkinBoundL0019:
    lemma_id: str = "L-0019"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    E0: float = 0.5
    nu: float = 0.1
    t: float = 0.02
    K2_star: int = 0
    S_t: float = 0.0
    Stokes_floor: float = 0.0
    envelope_cap: float = 0.0
    c0005_M: float = 61.23693461895651
    below_c0005: bool = False
    clay_implication: str = (
        "None. Stokes (linear) majorant on finite dealias Galerkin only; "
        "does not bound nonlinear NS; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0019(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    t: float = 0.02,
    c0005_M: float = 61.23693461895651,
) -> GalerkinBoundL0019:
    S, r_star = stokes_S(n, nu, t)
    floor = S * E0
    K2, _ = full_dealias_exact_stats(n)
    env = float(K2) * E0
    notes = (
        f"Stokes majorant S({t})={S:.6g} at |k|²={r_star}; "
        f"floor=S·E0={floor:.6g}. Envelope={env}. "
        f"floor < C-0005 M={c0005_M:.6g} < envelope: "
        f"{floor < c0005_M < env}. Nonlinear NS not controlled."
    )
    return GalerkinBoundL0019(
        n=n,
        E0=E0,
        nu=nu,
        t=t,
        K2_star=r_star,
        S_t=S,
        Stokes_floor=floor,
        envelope_cap=env,
        c0005_M=c0005_M,
        below_c0005=floor < c0005_M,
        notes=notes,
    )


def save_lemma_l0019(
    bound: GalerkinBoundL0019,
    path: str | Path = "conjectures/proved_restricted/L-0019.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
