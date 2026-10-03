"""
L-0026: Case-split / techo for C-0007 (gap below L-0024).

C-0007 asks Ω(T) ≤ M with M = mid(Stokes floor, L-0024 majorant) ≈ 41.284
on N≤24, E0=0.5, ν=0.1, T=0.02. The L-0024 comparison ODE is monotone in Ω0:
worst Ω(T) is at Ω0=E0 (≈41.743 > M); Ω(T) decreases to the Stokes floor as
Ω0 → K² E0.

Consequently there exists Ω★ ≈ 18.87 such that
  Ω0 ≥ Ω★  ⇒  Ω_L0024(T; Ω0) ≤ M.
The only obstruction to closing C-0007 with L-0024 is the low-Ω0 slab
  Ω0 ∈ [E0, Ω★).

Structural techo (N7 arithmetic + ODE probes):
- Replacing √(2 K² Ω) by a two-point upper stretch, clamping E_< ≤ min(E,δ/gap),
  or using E_K ≤ Ω/K² in the cross term, does **not** push the coupled majorant
  below M (endpoints stay ≥ ≈41.74 or jump to Young-like ≈45–50).
- min(defect, full Young) selects Young on the low slab and **worsens** Ω(T).
- Inflating the defect N-bound lowers the *numerical* ODE endpoint (faster
  overshoot → more viscous decay) but is **not** a valid way to claim a
  sharper majorant: the coupled comparison is not monotone in the size of N.

C-0007 remains open; needs a genuinely new low-Ω0 argument (or a larger M).

FINITE dealias Galerkin only. Not continuum. Not Clay.
Evidence: N7 (L-0024 reuse + case-split arithmetic).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0024_spectral_defect import (
    dealias_shell_stats,
    defect_ode_omega_T,
    worst_omega_T,
)
from ns_exploration.conjectures.l0025_multi_n_defect import c0007_open_target


def omega_star_for_M(
    M: float,
    E0: float,
    nu: float,
    T: float,
    stats: dict,
    tol: float = 1e-12,
    n_bisect: int = 48,
) -> tuple[float, float]:
    """
    Smallest Ω★ ∈ [E0, K² E0] with defect_ode_omega_T(E0, Ω★) ≤ M.
    Assumes Ω(T; Ω0) is nonincreasing in Ω0 (verified on the L-0024 grid).
    Returns (Ω★, Ω_L0024(T; Ω★)).
    """
    K2 = stats["K2"]
    lo = float(E0)
    hi = float(K2) * float(E0)
    if defect_ode_omega_T(E0, hi, nu, T, stats) > M + tol:
        raise RuntimeError("L-0024 at Ω0=K²E0 exceeds M; case-split impossible")
    if defect_ode_omega_T(E0, lo, nu, T, stats) <= M + tol:
        oT = defect_ode_omega_T(E0, lo, nu, T, stats)
        return lo, oT
    for _ in range(n_bisect):
        mid = 0.5 * (lo + hi)
        if defect_ode_omega_T(E0, mid, nu, T, stats) <= M + tol:
            hi = mid
        else:
            lo = mid
    oT = defect_ode_omega_T(E0, hi, nu, T, stats)
    return float(hi), float(oT)


def l0024_monotone_in_omega0(
    E0: float,
    nu: float,
    T: float,
    stats: dict,
    n_grid: int = 41,
) -> bool:
    xs = np.linspace(E0, stats["K2"] * E0, n_grid)
    vals = [defect_ode_omega_T(E0, float(o), nu, T, stats) for o in xs]
    return all(vals[i] + 1e-8 >= vals[i + 1] for i in range(len(vals) - 1))


@dataclass
class GalerkinBoundL0026:
    lemma_id: str = "L-0026"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    E0: float = 0.5
    nu: float = 0.1
    T: float = 0.02
    c0007_M: float = 0.0
    Stokes_floor: float = 0.0
    L0024_majorant: float = 0.0
    Omega_star: float = 0.0
    Omega_T_at_star: float = 0.0
    low_slab_L0024_worst: float = 0.0
    high_slab_closes: bool = False
    closes_c0007: bool = False
    monotone_in_Omega0: bool = False
    clay_implication: str = (
        "None. Case-split techo on L-0024 majorant only; C-0007 still open; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0026(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    n_grid: int = 61,
) -> GalerkinBoundL0026:
    gap = c0007_open_target(n=n, E0=E0, nu=nu, T=T)
    M = float(gap["proposed_bound_M"])
    stats = dealias_shell_stats(n)
    mono = l0024_monotone_in_omega0(E0, nu, T, stats, n_grid=min(n_grid, 41))
    om_star, omT_star = omega_star_for_M(M, E0, nu, T, stats)
    # Low-slab worst under L-0024 (= global worst when monotone)
    low_worst, _ = worst_omega_T(E0, nu, T, stats, n_grid=n_grid)
    # High slab: sample Ω0 ∈ [Ω★, K² E0]
    xs = np.linspace(om_star, stats["K2"] * E0, max(21, n_grid // 2))
    high_vals = [defect_ode_omega_T(E0, float(o), nu, T, stats) for o in xs]
    high_worst = float(max(high_vals))
    high_closes = high_worst <= M + 1e-9
    closes = low_worst <= M + 1e-9  # False for C-0007
    notes = (
        f"C-0007 M={M:.6f}; L-0024 majorant={gap['L0024_majorant']:.6f}; "
        f"Stokes floor={gap['Stokes_floor']:.6f}. "
        f"Monotone in Ω0: {mono}. Ω★={om_star:.6f} with Ω(T)={omT_star:.6f}≤M. "
        f"High slab [{om_star:.4f}, K²E0] worst={high_worst:.6f} closes={high_closes}. "
        f"Low slab still L-0024-worst={low_worst:.6f}>M; C-0007 open. "
        f"Techo: tightening/clamping N worsens or stays above M; "
        f"inflating N is not a valid sharper majorant."
    )
    return GalerkinBoundL0026(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        c0007_M=M,
        Stokes_floor=float(gap["Stokes_floor"]),
        L0024_majorant=float(gap["L0024_majorant"]),
        Omega_star=om_star,
        Omega_T_at_star=omT_star,
        low_slab_L0024_worst=float(low_worst),
        high_slab_closes=high_closes,
        closes_c0007=closes,
        monotone_in_Omega0=mono,
        notes=notes,
    )


def save_lemma_l0026(
    bound: GalerkinBoundL0026,
    path: str | Path = "conjectures/proved_restricted/L-0026.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
