"""
L-0002: Improved Galerkin enstrophy bound using ||∇u||_∞ ≤ √(3M)√(2Ω)
and viscous dissipation ||∇ω||² ≥ 2Ω (mean-zero modes |k|≥1).

ODE comparison:
  dΩ/dt ≤ -2ν Ω + C Ω^{3/2},  C = 2√2 √(3M),
reduces with z=√Ω, a=C/2=√2 √(3M) to
  z' = a z² - ν z,
with explicit solution (while denominator > 0):
  z(t) = ν / (a - (a - ν/z0) e^{ν t}),  z0 = √Ω(0) ≤ K√E0.

FINITE Galerkin ONLY. Not uniform in N. Not Clay. Evidence: N7.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0001_galerkin_bound import retained_mode_stats


@dataclass
class GalerkinBoundL0002:
    lemma_id: str = "L-0002"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    K: float = 0.0
    M: int = 0
    E0: float = 0.5
    nu: float = 0.1
    t: float = 0.02
    C_stretch: float = 0.0
    a: float = 0.0
    z0: float = 0.0
    Omega0_cap: float = 0.0
    Omega_t_cap: float | None = None
    finite_time_blowup_of_estimate: bool = False
    t_star_estimate: float | None = None
    clay_implication: str = (
        "None. M grows with resolution; estimate may even 'blow up' while true ODE stays smooth."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def stretch_constant(M: int) -> float:
    """C in dΩ/dt ≤ -2ν Ω + C Ω^{3/2}."""
    return 2.0 * math.sqrt(2.0) * math.sqrt(3.0 * M)


def explicit_bound_viscous(
    E0: float,
    t: float,
    K: float,
    M: int,
    nu: float = 0.1,
) -> GalerkinBoundL0002:
    """
    Ω(0) ≤ K² E0. Compare with z'=a z² - ν z, a = √2 √(3M).
    """
    if E0 <= 0 or M <= 0 or K < 0 or nu <= 0:
        raise ValueError("E0>0, M>0, K≥0, nu>0 required")
    C = stretch_constant(M)
    a = C / 2.0  # = √2 √(3M)
    Omega0 = (K**2) * E0
    z0 = math.sqrt(Omega0)
    # Denominator of z(t): a - (a - ν/z0) e^{ν t}
    # Blow-up of estimate if a z0 > ν and t large enough
    t_star = None
    blows = False
    if a * z0 > nu + 1e-15:
        t_star = (1.0 / nu) * math.log((a * z0) / (a * z0 - nu))
        if t >= t_star:
            blows = True

    if blows:
        omega_t = None
        notes = (
            f"Comparison ODE estimate blows at t*={t_star:.6e} ≤ t={t}; "
            "no finite upper bound from this differential inequality alone. "
            "Does NOT prove PDE/Galerkin blow-up."
        )
    else:
        denom = a - (a - nu / z0) * math.exp(nu * t)
        if denom <= 0:
            blows = True
            omega_t = None
            notes = "Denominator non-positive; estimate fails to close."
        else:
            z_t = nu / denom
            omega_t = z_t * z_t
            notes = (
                f"Viscous comparison: Ω(t)≤{omega_t:.6e} with C={C:.6g}, a={a:.6g}, "
                f"z0≤{z0:.6g}. Uses ||∇u||_∞≤√(3M)√(2Ω) and ||∇ω||²≥2Ω."
            )

    return GalerkinBoundL0002(
        K=K,
        M=M,
        E0=E0,
        nu=nu,
        t=t,
        C_stretch=C,
        a=a,
        z0=z0,
        Omega0_cap=Omega0,
        Omega_t_cap=omega_t,
        finite_time_blowup_of_estimate=blows,
        t_star_estimate=t_star,
        notes=notes,
    )


def bound_for_resolution_l0002(
    n: int,
    E0: float = 0.5,
    t: float = 0.02,
    nu: float = 0.1,
    dealias: bool = True,
) -> GalerkinBoundL0002:
    K, M = retained_mode_stats(n, dealias=dealias)
    b = explicit_bound_viscous(E0, t, K, M, nu=nu)
    b.notes += f" Grid N={n}, dealias={dealias}, K={K:.6g}, M={M}."
    return b


def save_lemma_l0002(bound: GalerkinBoundL0002, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
