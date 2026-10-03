"""
L-0010: Short-time Duhamel / bootstrap bound for shell-IC → full-grid Galerkin.

Under a bootstrap Ω(s) ≤ R on [0,T] and writing ε = (R - B)/K_full² for the
admissible high-mode energy (B = K_IC² E0), one has
  ||u||_∞ ≤ U_L + W √ε,
  ||u_H(t)||₂ ≤ ∫_0^t ||(u·∇)u||₂ ds ≤ T (U_L + W √ε) √(2 R),
hence E_H(t) ≤ T² (U_L + W √ε)² R.
The bootstrap closes when
  T² (U_L + W √ε)² R ≤ ε
and then Ω(t) ≤ R on [0,T].

Specializing to R = M of C-S-0002 gives an explicit horizon T_short on which
the C-S-0002 bound is *proved* for the finite Galerkin class (not the original
T=0.02 conjecture). Also reports the inviscid two-scale T_max from L-0009 and
a viscous RK4 diagnostic (N2 only).

FINITE Galerkin ONLY. Not continuum. Not Clay. Evidence: N7 (bootstrap), N2 (RK4).
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
from ns_exploration.conjectures.l0009_twoscale import (
    shell_mode_count,
    two_scale_cosh_cap,
)
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared, wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def kappa_H_min(n: int, k_ic: int, ic_kind: str) -> float | None:
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    k2 = k_squared(n)
    if ic_kind == "euclidean":
        low = mask & (np.sqrt(k2) <= float(k_ic) + 1e-12)
    elif ic_kind == "linf":
        linf = np.maximum(np.maximum(np.abs(kx), np.abs(ky)), np.abs(kz))
        low = mask & (linf <= float(k_ic) + 1e-12)
    else:
        raise ValueError(ic_kind)
    high = mask & ~low
    if not np.any(high):
        return None
    return float(np.min(k2[high]))


def duhamel_bootstrap_closes(
    B: float,
    K_full: float,
    U_L: float,
    W: float,
    R: float,
    T: float,
) -> tuple[bool, float, float]:
    """
    Returns (closes, epsilon, E_H_cap).
    Requires R > B and T ≥ 0.
    """
    if R <= B or T < 0:
        return False, 0.0, float("inf")
    eps = (R - B) / (K_full**2)
    if eps <= 0:
        return False, eps, float("inf")
    U = U_L + W * math.sqrt(eps)
    E_H_cap = (T**2) * (U**2) * R
    return bool(E_H_cap <= eps + 1e-15), eps, E_H_cap


def max_T_duhamel(
    B: float,
    K_full: float,
    U_L: float,
    W: float,
    R: float,
    T_hi: float = 0.05,
) -> float:
    """Largest T in [0, T_hi] for which the Duhamel bootstrap closes at level R."""
    if R <= B:
        return 0.0
    lo, hi = 0.0, T_hi
    # ensure hi fails or is cap
    if duhamel_bootstrap_closes(B, K_full, U_L, W, R, hi)[0]:
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if duhamel_bootstrap_closes(B, K_full, U_L, W, R, mid)[0]:
            lo = mid
        else:
            hi = mid
    return lo


def max_T_twoscale(
    B: float,
    K_full: float,
    U_L: float,
    W: float,
    T_hi: float = 0.05,
) -> float:
    lo, hi = 0.0, T_hi
    if two_scale_cosh_cap(B, K_full, U_L, W, hi)[0]:
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if two_scale_cosh_cap(B, K_full, U_L, W, mid)[0]:
            lo = mid
        else:
            hi = mid
    return lo


@dataclass
class GalerkinBoundL0010:
    lemma_id: str = "L-0010"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 16
    ic_kind: str = "linf"
    k_ic: int = 4
    K_ic_squared: int = 0
    K_full_squared: int = 0
    M_L: int = 0
    M_H: int = 0
    M_full: int = 0
    kappa_H: float | None = None
    E0: float = 0.5
    nu: float = 0.1
    B: float = 0.0
    U_L: float = 0.0
    W: float = 0.0
    # C-S-0002 short-time closure
    cs0002_M: float = 48.15928260103266
    T_duhamel_for_cs0002: float = 0.0
    proves_cs0002_on_short_horizon: bool = False
    T_target: float = 0.02
    duhamel_closes_at_T_target: bool = False
    # L-0009 horizon
    T_twoscale_max: float = 0.0
    twoscale_cap_at_Tmax: float | None = None
    clay_implication: str = (
        "None. Proves a SHORT-TIME finite-Galerkin bound for the C-S IC class; "
        "does not prove C-S-0002 at T=0.02 and has no Clay implication."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0010(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    nu: float = 0.1,
    T_target: float = 0.02,
    cs0002_M: float = 48.15928260103266,
) -> GalerkinBoundL0010:
    if ic_kind == "linf":
        K_ic2 = linf_shell_K2(k_ic)
    else:
        K_ic2 = euclidean_shell_K2(n, k_ic)

    K2_full, M_full = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = shell_mode_count(n, k_ic, ic_kind)
    M_H = M_full - M_L
    B = float(K_ic2) * E0
    U_L = math.sqrt(3.0 * max(M_L, 1)) * math.sqrt(2.0 * E0)
    W = math.sqrt(3.0 * M_H) * math.sqrt(2.0) if M_H > 0 else 0.0
    kappa = kappa_H_min(n, k_ic, ic_kind)

    T_duh = max_T_duhamel(B, K_full, U_L, W, cs0002_M)
    closes_target = duhamel_bootstrap_closes(B, K_full, U_L, W, cs0002_M, T_target)[0]
    T_ts = max_T_twoscale(B, K_full, U_L, W)
    cap_ts = None
    if T_ts > 0:
        ok, cap_ts, _, _ = two_scale_cosh_cap(B, K_full, U_L, W, T_ts * 0.999)
        if not ok:
            cap_ts = None

    notes = (
        f"IC {ic_kind} k={k_ic} N={n}: B={B:.6g}, U_L={U_L:.6g}, W={W:.6g}, "
        f"kappa_H={kappa}. Duhamel T* for R={cs0002_M:.4g}: {T_duh:.6g} "
        f"(proves C-S bound on [0,T*]). At T={T_target}: closes={closes_target}. "
        f"Two-scale T_max={T_ts:.6g}."
    )
    return GalerkinBoundL0010(
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        M_full=M_full,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        B=B,
        U_L=U_L,
        W=W,
        cs0002_M=cs0002_M,
        T_duhamel_for_cs0002=T_duh,
        proves_cs0002_on_short_horizon=bool(T_duh > 0),
        T_target=T_target,
        duhamel_closes_at_T_target=closes_target,
        T_twoscale_max=T_ts,
        twoscale_cap_at_Tmax=cap_ts,
        notes=notes,
    )


def save_lemma_l0010(bound: GalerkinBoundL0010, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
