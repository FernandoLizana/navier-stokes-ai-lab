"""
L-0003: Uniform-in-time Galerkin enstrophy bound via enstrophy–energy Cauchy.

Key estimate (mean-zero Fourier–Galerkin on T³):
  Σ |k|⁴ |û_k|² ≥ (Σ |k|² |û_k|²)² / (Σ |û_k|²) = 2 Ω² / E,
so ||∇ω||₂² ≥ 2 Ω² / E ≥ 2 Ω² / E0  (using E(t)≤E0).

With ||∇u||_∞ ≤ √(3M) √(2Ω) as in L-0002,
  dΩ/dt ≤ -ν (2 Ω² / E0) + C Ω^{3/2},   C = 2√2 √(3M).

Comparison ODE for z=√Ω:
  z' = a z² - b z³,   a = C/2 = √2 √(3M),   b = ν/E0,
has globally attracting equilibrium z_eq = a/b (if a,b>0), hence
  Ω(t) ≤ max( Ω(0), Ω_eq ) ≤ max( K² E0,  6 M E0² / ν² )
for all t ≥ 0.

FINITE Galerkin ONLY. Not uniform as N→∞ (M→∞). Not Clay. Evidence: N7.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0001_galerkin_bound import retained_mode_stats
from ns_exploration.conjectures.l0002_viscous_galerkin import stretch_constant


@dataclass
class GalerkinBoundL0003:
    lemma_id: str = "L-0003"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    K: float = 0.0
    M: int = 0
    E0: float = 0.5
    nu: float = 0.1
    Omega0_cap: float = 0.0
    Omega_eq: float = 0.0
    Omega_uniform_cap: float = 0.0
    closes: bool = True
    clay_implication: str = (
        "None. Uniform in time for fixed N only; Ω_eq ∼ M/ν² → ∞ as resolution → ∞."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def explicit_uniform_bound(
    E0: float,
    K: float,
    M: int,
    nu: float = 0.1,
) -> GalerkinBoundL0003:
    if E0 <= 0 or M <= 0 or nu <= 0 or K < 0:
        raise ValueError("require E0>0, M>0, nu>0, K≥0")
    C = stretch_constant(M)
    a = C / 2.0
    b = nu / E0
    Omega0 = (K**2) * E0
    # z_eq = a/b ⇒ Ω_eq = (a/b)² = a² E0² / ν²
    # a² = 2 * 3M = 6M ⇒ Ω_eq = 6 M E0² / ν²
    Omega_eq = 6.0 * M * (E0**2) / (nu**2)
    Omega_cap = max(Omega0, Omega_eq)
    return GalerkinBoundL0003(
        K=K,
        M=M,
        E0=E0,
        nu=nu,
        Omega0_cap=Omega0,
        Omega_eq=Omega_eq,
        Omega_uniform_cap=Omega_cap,
        closes=True,
        notes=(
            f"Uniform-in-time Galerkin bound: Ω(t)≤max(K²E0, 6M E0²/ν²)="
            f"{Omega_cap:.6e} (Ω_eq={Omega_eq:.6e}, Ω0≤{Omega0:.6e}). "
            f"C={C:.6g}, a={a:.6g}, b={b:.6g}."
        ),
    )


def bound_for_resolution_l0003(
    n: int,
    E0: float = 0.5,
    nu: float = 0.1,
    dealias: bool = True,
) -> GalerkinBoundL0003:
    K, M = retained_mode_stats(n, dealias=dealias)
    b = explicit_uniform_bound(E0, K, M, nu=nu)
    b.notes += f" Grid N={n}, dealias={dealias}, K={K:.6g}, M={M}."
    return b


def save_lemma_l0003(bound: GalerkinBoundL0003, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
