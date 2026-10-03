"""
L-0001: Explicit (crude) a priori enstrophy bound for Fourier-Galerkin NS on T³.

Hypotheses (finite-dimensional ONLY):
  - Divergence-free velocity expanded in finitely many Fourier modes
  - Every active wavevector satisfies |k| ≤ K (Euclidean)
  - At most M complex mode slots counted in the ℓ¹/Cauchy–Schwarz embedding
    (safe: M = (# retained lattice points))
  - Kinetic energy E(t) ≤ E0 for all t (true for f=0 since dE/dt = -ν||∇u||² ≤ 0)
  - Domain T³ = [0,2π)³ with NS-MRL Fourier conventions

Conclusion:
  Ω(t) ≤ K² E0 * exp(2 K √(3M) √(2 E0) t)
where Ω = (1/2)||ω||₂² (enstrophy) and M = # retained spatial lattice points
(vector field: Cauchy–Schwarz over 3M complex coefficients).

Proof sketch (hand):
  1) Ω = (1/2) Σ |k|²|û_k|² ≤ (1/2) K² Σ|û_k|² = K² E.
  2) dΩ/dt = -ν||∇ω||₂² + ∫ ω·(ω·∇)u ≤ ∫ ω·(ω·∇)u
     (drop dissipation for an upper bound).
  3) |∫ ω·(ω·∇)u| ≤ ||ω||₂² ||∇u||_∞ = 2 Ω ||∇u||_∞.
  4) ||∇u||_∞ ≤ Σ_{i,k} |k| |û_k^i| ≤ K Σ_{i,k} |û_k^i|
     ≤ K √(3M) (Σ_{i,k}|û_k^i|²)^{1/2} = K √(3M) √(2E) ≤ K √(3M) √(2 E0).
  5) Hence dΩ/dt ≤ 2 K √(3M) √(2 E0) Ω, so Ω(t) ≤ Ω(0) e^{α t}
     with α = 2 K √(3M) √(2 E0) and Ω(0) ≤ K² E0.

This is N7 for the *finite Galerkin ODE* under the stated hypotheses.
It is NOT uniform in N→∞ (K,M →∞), NOT a continuum bound, NOT Clay.

Evidence level: N7 (elementary estimates). Numerical checks are N2 confirmations only.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared, wave_number_grids


@dataclass
class GalerkinEnstrophyBound:
    lemma_id: str = "L-0001"
    route: str = "B"
    status: str = "proved_restricted"  # finite Galerkin only
    evidence_level: str = "N7"
    K: float = 0.0
    M: int = 0
    E0: float = 0.5
    t: float = 0.02
    alpha: float = 0.0
    Omega0_cap: float = 0.0
    Omega_t_cap: float = 0.0
    clay_implication: str = (
        "None. Constants blow up as K,M→∞; no continuum regularity."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def retained_mode_stats(n: int, dealias: bool = True) -> tuple[float, int]:
    """
    Return (K_eucl_max, M_count) for the computational grid.
    M counts lattice points retained by the dealias mask (all k, including 0).
    """
    kx, ky, kz = wave_number_grids(n)
    k2 = k_squared(n)
    if dealias:
        mask = dealias_mask(n)
    else:
        mask = np.ones_like(k2, dtype=bool)
    k2_ret = k2[mask]
    K = float(np.sqrt(np.max(k2_ret))) if k2_ret.size else 0.0
    M = int(np.count_nonzero(mask))
    return K, M


def explicit_enstrophy_bound(E0: float, t: float, K: float, M: int) -> GalerkinEnstrophyBound:
    alpha = 2.0 * K * math.sqrt(3.0 * M) * math.sqrt(2.0 * E0)
    Omega0 = (K**2) * E0
    Omega_t = Omega0 * math.exp(alpha * t)
    return GalerkinEnstrophyBound(
        K=K,
        M=M,
        E0=E0,
        t=t,
        alpha=alpha,
        Omega0_cap=Omega0,
        Omega_t_cap=Omega_t,
        notes=(
            "Crude exponential bound dropping viscosity; valid for Galerkin with "
            f"|k|≤K and M spatial modes (√(3M) vector CS). For E0={E0}, t={t}: "
            f"Ω(t)≤{Omega_t:.6e}."
        ),
    )


def bound_for_resolution(n: int, E0: float = 0.5, t: float = 0.02, dealias: bool = True) -> GalerkinEnstrophyBound:
    K, M = retained_mode_stats(n, dealias=dealias)
    b = explicit_enstrophy_bound(E0, t, K, M)
    b.notes += f" Grid N={n}, dealias={dealias}, K={K:.6g}, M={M}."
    return b


def save_lemma(bound: GalerkinEnstrophyBound, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
