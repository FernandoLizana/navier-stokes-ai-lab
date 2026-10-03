"""
L-0008: High-mode cascade bound for shell ICs on the full dealias grid.

Key improvement over L-0007: uses E_H(0) = 0 (high modes start empty), not only
the shell enstrophy ceiling.

Split the dealias lattice into
  L = {k : |k| ≤ K_IC} (or |k|_∞ ≤ k_inf),   H = dealias \\ L.
Fourier orthogonality gives E = E_L + E_H, Ω = Ω_L + Ω_H, and for all t
  Ω_L(t) ≤ K_IC² E_L(t) ≤ K_IC² E0 =: B,
  Ω_H(t) ≤ K_full² E_H(t).
Hence Ω(t) ≤ B + K_full² E_H(t).

Energy of high modes (drop viscosity for an upper bound):
  dE_H/dt ≤ ||u_H||₂ ||(u·∇)u||₂ ≤ √(2 E_H) · ||u||_∞ · √(2Ω).
With ||u||_∞ ≤ √(3 M_full) √(2 E) ≤ √(3 M_full) √(2 E0),
  d√E_H / dt ≤ √2 √(3 M_full) √E0 · √Ω =: A √Ω.
Using Ω ≤ B + K_full² E_H and the comparison ODE
  z' = A √(B + K_full² z²),  z(0)=0,  z=√E_H,
one obtains z(t) = √(B)/K_full · sinh(A K_full t), therefore
  Ω(t) ≤ B · cosh²(A K_full t),   A = √2 √(3 M_full) √E0.

FINITE Galerkin ONLY. Not continuum. Not Clay. Evidence: N7.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0003_uniform_galerkin import explicit_uniform_bound
from ns_exploration.conjectures.l0007_shell_ic_fullgrid import (
    euclidean_shell_K2,
    lemma_l0007,
    linf_shell_K2,
)
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class GalerkinBoundL0008:
    lemma_id: str = "L-0008"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 12
    ic_kind: str = "euclidean"
    k_ic: int = 2
    K_ic_squared: int = 0
    K_full_squared: int = 0
    M_full: int = 0
    E0: float = 0.5
    nu: float = 0.1
    t: float = 0.02
    Omega0_shell: float = 0.0
    A_cascade: float = 0.0
    alpha: float = 0.0  # A * K_full
    cosh_cap: float = 0.0
    L0007_cap: float = 0.0
    L0003_cap: float = 0.0
    best_cap: float = 0.0
    which_best: str = ""
    proves_cs0002: bool = False
    cs0002_M: float | None = None
    clay_implication: str = (
        "None. Finite Galerkin; uses E_H(0)=0 but M_full stretch; not continuum NS."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def cascade_cosh_bound(
    Omega0_shell: float,
    M_full: int,
    K_full: float,
    E0: float,
    t: float,
) -> tuple[float, float, float]:
    """
    Returns (A, alpha, Omega_cap) with
      A = √2 √(3M) √E0,  alpha = A K_full,  Ω ≤ Ω0 cosh²(alpha t).
    """
    A = math.sqrt(2.0) * math.sqrt(3.0 * M_full) * math.sqrt(E0)
    alpha = A * K_full
    # Use cosh via exp for stability: cosh(x) = (e^x + e^{-x})/2
    x = alpha * t
    if x > 50:
        # Overflow guard: bound is huge; report +inf-like large float
        cosh_x = 0.5 * math.exp(x)
    else:
        cosh_x = math.cosh(x)
    return A, alpha, Omega0_shell * (cosh_x**2)


def lemma_l0008(
    n: int = 12,
    k_ic: int = 2,
    ic_kind: str = "euclidean",
    E0: float = 0.5,
    t: float = 0.02,
    nu: float = 0.1,
    cs0002_M: float | None = 48.15928260103266,
) -> GalerkinBoundL0008:
    if ic_kind == "linf":
        K_ic2 = linf_shell_K2(k_ic)
    elif ic_kind == "euclidean":
        K_ic2 = euclidean_shell_K2(n, k_ic)
    else:
        raise ValueError(ic_kind)

    K2_full, M_full = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    B = float(K_ic2) * E0

    A, alpha, cosh_cap = cascade_cosh_bound(B, M_full, K_full, E0, t)
    b7 = lemma_l0007(
        n=n, k_ic=k_ic, ic_kind=ic_kind, E0=E0, t=t, nu=nu, cs0002_M=cs0002_M
    )
    b3 = explicit_uniform_bound(E0, K_full, M_full, nu=nu)

    candidates = [
        ("L-0008-cosh", cosh_cap),
        ("L-0007", b7.best_cap),
        ("L-0003-full", b3.Omega_uniform_cap),
    ]
    which, best = min(candidates, key=lambda x: x[1])
    proves = bool(cs0002_M is not None and best <= float(cs0002_M))

    notes = (
        f"IC {ic_kind} k={k_ic} N={n}: B=K_IC^2 E0={B:.6g}, M={M_full}, "
        f"K_full={K_full:.6g}, A={A:.6g}, alpha={alpha:.6g}, "
        f"cosh_cap={cosh_cap:.6e}. Best={which} -> {best:.6e}. "
        f"Proves C-S-0002: {proves}."
    )
    return GalerkinBoundL0008(
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_full=M_full,
        E0=E0,
        nu=nu,
        t=t,
        Omega0_shell=B,
        A_cascade=A,
        alpha=alpha,
        cosh_cap=cosh_cap,
        L0007_cap=b7.best_cap,
        L0003_cap=b3.Omega_uniform_cap,
        best_cap=best,
        which_best=which,
        proves_cs0002=proves,
        cs0002_M=cs0002_M,
        notes=notes,
    )


def save_lemma_l0008(bound: GalerkinBoundL0008, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
