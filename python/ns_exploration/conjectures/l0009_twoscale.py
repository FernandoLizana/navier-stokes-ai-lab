"""
L-0009: Two-scale ||u||_∞ embedding for shell IC → full dealias grid.

Split dealias modes into L (IC shell) and H (complement):
  ||u||_∞ ≤ √(3 M_L) √(2 E_L) + √(3 M_H) √(2 E_H)
           ≤ U_L + W √E_H,
  U_L = √(3 M_L) √(2 E0),  W = √(3 M_H) √2  (E_L ≤ E0).

With Ω ≤ B + K_full² E_H, B = K_IC² E0, z = √E_H, z(0)=0:
  z' ≤ (U_L + W z) √(B + K_full² z²).

Hyperbolic substitution z = (√B / K_full) sinh u yields
  t = (1/K_full) ∫_0^u ds / (U_L + λ sinh s),   λ = W √B / K_full,
and the explicit antiderivative
  F(u) = (1/√(U_L²+λ²)) log |(τ - λ/U_L + σ)/(τ - λ/U_L - σ)|,
  τ = tanh(u/2), σ = √(U_L²+λ²)/U_L,
so Ω(t) = B cosh²(u(t)). If the integral saturates before t (F(∞)-F(0) < K t),
the comparison ODE has no global solution on [0,t]; fall back to L-0008/L-0007/L-0003.

FINITE Galerkin ONLY. Not Clay. Evidence: N7.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import (
    euclidean_shell_K2,
    linf_shell_K2,
)
from ns_exploration.conjectures.l0008_cascade import lemma_l0008
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared, wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def shell_mode_count(n: int, k_ic: int, ic_kind: str) -> int:
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    if ic_kind == "euclidean":
        kr = np.sqrt(k_squared(n))
        sm = mask & (kr <= float(k_ic) + 1e-12)
    elif ic_kind == "linf":
        linf = np.maximum(np.maximum(np.abs(kx), np.abs(ky)), np.abs(kz))
        sm = mask & (linf <= float(k_ic) + 1e-12)
    else:
        raise ValueError(ic_kind)
    return int(sm.sum())


def _F_of_u(u: float, U_L: float, lam: float) -> float:
    """Antiderivative F(u) for ∫ du/(U_L + λ sinh u)."""
    if U_L <= 0:
        raise ValueError("U_L must be positive")
    s = math.sqrt(U_L * U_L + lam * lam)
    th = math.tanh(0.5 * u)
    num = th - lam / U_L + s / U_L
    den = th - lam / U_L - s / U_L
    if abs(den) < 1e-300:
        return float("inf")
    return math.log(abs(num / den)) / s


def two_scale_cosh_cap(
    B: float,
    K_full: float,
    U_L: float,
    W: float,
    t: float,
) -> tuple[bool, float | None, float | None, float]:
    """
    Returns (closes, Omega_cap, u_star, F_inf_minus_F0).
    """
    if t < 0:
        raise ValueError("t>=0")
    if t == 0.0:
        return True, B, 0.0, float("inf")
    if W <= 0.0:
        # Reduces to L-0008: u = U_L K_full t, Ω = B cosh²(u)
        u = U_L * K_full * t
        return True, B * (math.cosh(u) ** 2), u, float("inf")

    lam = W * math.sqrt(B) / K_full
    F0 = _F_of_u(0.0, U_L, lam)
    # F(∞): tanh→1 (use same abs-log formula as _F_of_u)
    s = math.sqrt(U_L * U_L + lam * lam)
    num_inf = 1.0 - lam / U_L + s / U_L
    den_inf = 1.0 - lam / U_L - s / U_L
    if abs(den_inf) < 1e-300 or abs(num_inf) < 1e-300:
        # integral diverges ⇒ comparison always closes for all t
        F_inf = float("inf")
    else:
        F_inf = math.log(abs(num_inf / den_inf)) / s
    span = F_inf - F0
    need = K_full * t
    if math.isfinite(span) and need >= span - 1e-15:
        return False, None, None, span

    target = need + F0
    lo, hi = 0.0, 1.0
    # expand hi until F(hi) >= target
    for _ in range(60):
        if _F_of_u(hi, U_L, lam) >= target:
            break
        hi *= 2.0
    else:
        return False, None, None, span

    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if _F_of_u(mid, U_L, lam) < target:
            lo = mid
        else:
            hi = mid
    u_star = hi
    return True, B * (math.cosh(u_star) ** 2), u_star, span


@dataclass
class GalerkinBoundL0009:
    lemma_id: str = "L-0009"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 12
    ic_kind: str = "euclidean"
    k_ic: int = 2
    K_ic_squared: int = 0
    K_full_squared: int = 0
    M_L: int = 0
    M_H: int = 0
    M_full: int = 0
    E0: float = 0.5
    t: float = 0.02
    nu: float = 0.1
    Omega0_shell: float = 0.0
    U_L: float = 0.0
    W: float = 0.0
    two_scale_closes: bool = False
    two_scale_cap: float | None = None
    u_star: float | None = None
    L0008_cap: float = 0.0
    best_cap: float = 0.0
    which_best: str = ""
    proves_cs0002: bool = False
    cs0002_M: float | None = None
    clay_implication: str = (
        "None. Finite Galerkin two-scale embedding; not continuum NS."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0009(
    n: int = 12,
    k_ic: int = 2,
    ic_kind: str = "euclidean",
    E0: float = 0.5,
    t: float = 0.02,
    nu: float = 0.1,
    cs0002_M: float | None = 48.15928260103266,
) -> GalerkinBoundL0009:
    if ic_kind == "linf":
        K_ic2 = linf_shell_K2(k_ic)
    else:
        K_ic2 = euclidean_shell_K2(n, k_ic)

    K2_full, M_full = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = shell_mode_count(n, k_ic, ic_kind)
    M_H = M_full - M_L
    B = float(K_ic2) * E0
    U_L = math.sqrt(3.0 * M_L) * math.sqrt(2.0 * E0)
    W = math.sqrt(3.0 * M_H) * math.sqrt(2.0) if M_H > 0 else 0.0

    closes, cap9, u_star, _span = two_scale_cosh_cap(B, K_full, U_L, W, t)
    b8 = lemma_l0008(
        n=n, k_ic=k_ic, ic_kind=ic_kind, E0=E0, t=t, nu=nu, cs0002_M=cs0002_M
    )

    candidates: list[tuple[str, float]] = [
        ("L-0008", b8.best_cap),
        ("L-0007", b8.L0007_cap),
        ("L-0003-full", b8.L0003_cap),
    ]
    if closes and cap9 is not None:
        candidates.append(("L-0009-two-scale", cap9))

    which, best = min(candidates, key=lambda x: x[1])
    proves = bool(cs0002_M is not None and best <= float(cs0002_M))

    notes = (
        f"IC {ic_kind} k={k_ic} N={n}: M_L={M_L}, M_H={M_H}, B={B:.6g}, "
        f"U_L={U_L:.6g}, W={W:.6g}, two_scale_closes={closes}, "
        f"cap9={cap9}. Best={which} -> {best:.6e}. Proves C-S-0002: {proves}."
    )
    return GalerkinBoundL0009(
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        M_full=M_full,
        E0=E0,
        t=t,
        nu=nu,
        Omega0_shell=B,
        U_L=U_L,
        W=W,
        two_scale_closes=closes,
        two_scale_cap=cap9,
        u_star=u_star,
        L0008_cap=b8.best_cap,
        best_cap=best,
        which_best=which,
        proves_cs0002=proves,
        cs0002_M=cs0002_M,
        notes=notes,
    )


def save_lemma_l0009(bound: GalerkinBoundL0009, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
