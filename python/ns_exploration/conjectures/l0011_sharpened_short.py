"""
L-0011: Sharper short-time Duhamel bounds for the C-S class.

Improvements over L-0010:
  1) Divergence-free Fourier embedding: ||u||_∞ ≤ √(2 M) √(2 E)
     (2 complex dof per mode, not 3).
  2) Viscous Duhamel factor τ(T) = (1 - e^{-ν κ_H T}) / (ν κ_H)
     instead of the crude length T.
  3) L×L constant-force majorant: if one bounds only the shell–shell
     contribution to high-mode production by U_L √(2 B), then
       E_H(t) ≤ τ(T)² U_L² B,
       Ω(t) ≤ B (1 + (K_full U_L τ(T))²).
     This is recorded as a *conditional* majorant (N6) unless the
     H-cross terms are absorbed; the fully rigorous package is (1)+(2)
     inside the L-0010-style bootstrap (N7).

Outputs:
  - T*_rigorous for R = C-S-0002 M  (div-free + viscous bootstrap)
  - T*_LxL for the same R (conditional majorant)
  - Ω_LxL(T=0.02) as a conditional a priori number for the full horizon

FINITE Galerkin ONLY. Not Clay.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2
from ns_exploration.conjectures.l0009_twoscale import shell_mode_count
from ns_exploration.conjectures.l0010_duhamel_short import kappa_H_min
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def embedding_U(M: int, E0: float, div_free: bool = True) -> float:
    """||u||_∞ ≤ U with kinetic energy ≤ E0 on M lattice modes."""
    if M <= 0:
        return 0.0
    dof = 2.0 if div_free else 3.0
    return math.sqrt(dof * M) * math.sqrt(2.0 * E0)


def embedding_W(M_H: int, div_free: bool = True) -> float:
    """||u_H||_∞ ≤ W √E_H."""
    if M_H <= 0:
        return 0.0
    dof = 2.0 if div_free else 3.0
    # √(dof M_H) √(2 E_H) = √(2 dof M_H) √E_H
    return math.sqrt(2.0 * dof * M_H)


def tau_viscous(T: float, nu: float, kappa_H: float) -> float:
    if T <= 0:
        return 0.0
    a = nu * kappa_H
    if a <= 0:
        return T
    return (1.0 - math.exp(-a * T)) / a


def bootstrap_closes(
    B: float,
    K_full: float,
    U_L: float,
    W: float,
    R: float,
    T: float,
    nu: float,
    kappa_H: float,
) -> tuple[bool, float, float, float]:
    """Returns (closes, eps, E_H_cap, tau)."""
    if R <= B or T < 0:
        return False, 0.0, float("inf"), 0.0
    eps = (R - B) / (K_full**2)
    U = U_L + W * math.sqrt(eps)
    tauf = tau_viscous(T, nu, kappa_H)
    E_H_cap = (tauf**2) * (U**2) * R
    return bool(E_H_cap <= eps + 1e-15), eps, E_H_cap, tauf


def max_T_bootstrap(
    B: float,
    K_full: float,
    U_L: float,
    W: float,
    R: float,
    nu: float,
    kappa_H: float,
    T_hi: float = 0.05,
) -> float:
    lo, hi = 0.0, T_hi
    if bootstrap_closes(B, K_full, U_L, W, R, hi, nu, kappa_H)[0]:
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if bootstrap_closes(B, K_full, U_L, W, R, mid, nu, kappa_H)[0]:
            lo = mid
        else:
            hi = mid
    return lo


def omega_LxL_majorant(
    B: float,
    K_full: float,
    U_L: float,
    T: float,
    nu: float,
    kappa_H: float,
) -> float:
    """Conditional: Ω ≤ B (1 + (K U_L τ)^2) from L×L-only forcing."""
    tauf = tau_viscous(T, nu, kappa_H)
    return B * (1.0 + (K_full * U_L * tauf) ** 2)


def max_T_LxL(
    B: float,
    K_full: float,
    U_L: float,
    R: float,
    nu: float,
    kappa_H: float,
    T_hi: float = 0.05,
) -> float:
    lo, hi = 0.0, T_hi
    if omega_LxL_majorant(B, K_full, U_L, hi, nu, kappa_H) <= R:
        return hi
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if omega_LxL_majorant(B, K_full, U_L, mid, nu, kappa_H) <= R:
            lo = mid
        else:
            hi = mid
    return lo


@dataclass
class GalerkinBoundL0011:
    lemma_id: str = "L-0011"
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
    kappa_H: float = 0.0
    E0: float = 0.5
    nu: float = 0.1
    B: float = 0.0
    U_L: float = 0.0
    W: float = 0.0
    cs0002_M: float = 48.15928260103266
    T_target: float = 0.02
    # Rigorous bootstrap (div-free + viscous)
    T_star_rigorous: float = 0.0
    improvement_vs_L0010: float = 0.0
    # Conditional L×L majorant
    T_star_LxL_conditional: float = 0.0
    Omega_LxL_at_T_target: float = 0.0
    LxL_evidence: str = "N6"
    clay_implication: str = (
        "None. Finite Galerkin short-time bounds only; C-S-0002 at T=0.02 still open."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0011(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    nu: float = 0.1,
    T_target: float = 0.02,
    cs0002_M: float = 48.15928260103266,
    T_L0010_ref: float = 0.0010112385074348496,
) -> GalerkinBoundL0011:
    if ic_kind != "linf":
        raise ValueError("L-0011 packaged for linf C-S class")

    K_ic2 = linf_shell_K2(k_ic)
    K2_full, M_full = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = shell_mode_count(n, k_ic, ic_kind)
    M_H = M_full - M_L
    kappa = kappa_H_min(n, k_ic, ic_kind) or 0.0
    B = float(K_ic2) * E0
    U_L = embedding_U(M_L, E0, div_free=True)
    W = embedding_W(M_H, div_free=True)

    T_rig = max_T_bootstrap(B, K_full, U_L, W, cs0002_M, nu, kappa)
    T_lxl = max_T_LxL(B, K_full, U_L, cs0002_M, nu, kappa)
    om_tgt = omega_LxL_majorant(B, K_full, U_L, T_target, nu, kappa)
    improv = T_rig / T_L0010_ref if T_L0010_ref > 0 else float("inf")

    notes = (
        f"div-free+viscous bootstrap T*={T_rig:.8f} for R={cs0002_M:.4g} "
        f"({improv:.3g}x vs L-0010). Conditional LxL T*={T_lxl:.8f}; "
        f"LxL Omega(T={T_target})={om_tgt:.6g} (N6, ignores H-bilinear remainder)."
    )
    return GalerkinBoundL0011(
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        B=B,
        U_L=U_L,
        W=W,
        cs0002_M=cs0002_M,
        T_target=T_target,
        T_star_rigorous=T_rig,
        improvement_vs_L0010=improv,
        T_star_LxL_conditional=T_lxl,
        Omega_LxL_at_T_target=om_tgt,
        notes=notes,
    )


def save_lemma_l0011(bound: GalerkinBoundL0011, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
