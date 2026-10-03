"""
L-0006: Improved short-time shell Galerkin bound.

Stretching via Cauchy–Schwarz on coefficients (no extra K factor):
  ||∇u||_∞ ≤ Σ_{i,k} |k| |û_k^i| ≤ √(Σ |k|²|û|²) √(3M) = √(2Ω) √(3M),
hence
  dΩ/dt ≤ 2√2 √(3M) Ω^{3/2} =: 2 a Ω^{3/2},  a = √2 √(3M).

Comparison (inviscid upper): z=√Ω satisfies z' ≤ a z², so while t < 1/(a z0),
  Ω(t) ≤ Ω0 / (1 - a √Ω0 t)²,   Ω0 ≤ K² E0.

When that denominator would vanish on [0,T], fall back to L-0001 exponential
(which keeps a factor K√E0 instead of √Ω and can remain finite longer).

Also report L-0003 uniform ceiling. Best finite bound among those that close.

FINITE shell Galerkin ONLY. Not Clay. Evidence: N7.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0001_galerkin_bound import explicit_enstrophy_bound
from ns_exploration.conjectures.l0003_uniform_galerkin import explicit_uniform_bound
from ns_exploration.conjectures.l0004_spectral_support import shell_a_priori_bounds


@dataclass
class GalerkinBoundL0006:
    lemma_id: str = "L-0006"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 12
    k_shell: float = 3.0
    K: float = 0.0
    M: int = 0
    E0: float = 0.5
    nu: float = 0.1
    t: float = 0.02
    a_stretch: float = 0.0
    Omega0_cap: float = 0.0
    algebraic_closes: bool = False
    algebraic_cap: float | None = None
    t_star_algebraic: float | None = None
    L0001_cap: float = 0.0
    L0003_cap: float = 0.0
    best_cap: float = 0.0
    which_best: str = ""
    clay_implication: str = "None. Shell Galerkin only; constants grow with M."
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def algebraic_short_time_bound(E0: float, t: float, K: float, M: int) -> tuple[bool, float | None, float, float]:
    """
    Returns (closes, Omega_t, a, t_star).
    Ω(t) ≤ Ω0/(1-a√Ω0 t)² with Ω0=K²E0, a=√2√(3M).
    """
    a = math.sqrt(2.0) * math.sqrt(3.0 * M)
    Omega0 = (K**2) * E0
    z0 = math.sqrt(Omega0)
    t_star = 1.0 / (a * z0) if a * z0 > 0 else float("inf")
    if t >= t_star:
        return False, None, a, t_star
    denom = 1.0 - a * z0 * t
    omega_t = Omega0 / (denom**2)
    return True, omega_t, a, t_star


def lemma_l0006(
    n: int = 12,
    k_shell: float = 3.0,
    E0: float = 0.5,
    t: float = 0.02,
    nu: float = 0.1,
) -> GalerkinBoundL0006:
    shell = shell_a_priori_bounds(n, int(k_shell), E0=E0, t=t, nu=nu)
    K, M = shell.K, shell.M
    closes, alg, a, t_star = algebraic_short_time_bound(E0, t, K, M)
    b1 = explicit_enstrophy_bound(E0, t, K, M)
    b3 = explicit_uniform_bound(E0, K, M, nu=nu)

    candidates: list[tuple[str, float]] = [("L-0001", b1.Omega_t_cap), ("L-0003", b3.Omega_uniform_cap)]
    if closes and alg is not None:
        candidates.append(("L-0006-algebraic", alg))

    which, best = min(candidates, key=lambda x: x[1])
    notes = (
        f"Shell |k|≤{k_shell} on N={n}: K={K:.6g}, M={M}, a={a:.6g}, t*={t_star:.6g}. "
        f"Algebraic closes={closes}. Best={which} → {best:.6e}."
    )
    return GalerkinBoundL0006(
        n=n,
        k_shell=k_shell,
        K=K,
        M=M,
        E0=E0,
        nu=nu,
        t=t,
        a_stretch=a,
        Omega0_cap=(K**2) * E0,
        algebraic_closes=closes,
        algebraic_cap=alg,
        t_star_algebraic=t_star,
        L0001_cap=b1.Omega_t_cap,
        L0003_cap=b3.Omega_uniform_cap,
        best_cap=best,
        which_best=which,
        notes=notes,
    )


def save_l0006(b: GalerkinBoundL0006, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(b.as_dict(), indent=2), encoding="utf-8")
