"""
L-0012: High-mode energy ODE using Navier–Stokes cancellation identities.

For the Leray-projected Galerkin system, with L = IC shell and H = complement,
  dE_H/dt = -ν ||∇u_H||₂² - ⟨u_H, (u·∇)u⟩.
The identities
  ⟨u_H, (u_L·∇)u_H⟩ = 0,   ⟨u_H, (u_H·∇)u_H⟩ = 0
(for divergence-free fields) leave only
  ⟨u_H, (u_L·∇)u_L⟩ + ⟨u_H, (u_H·∇)u_L⟩.
Hence, with Ω_L ≤ B = K_IC² E0, ||u_L||_∞ ≤ U_L, ||u_H||_∞ ≤ W √E_H,
  κ_H = min_{k∈H} |k|²,
  dE_H/dt ≤ -2 ν κ_H E_H + 2 U_L √B √E_H + 2 W √B E_H.
For z = √E_H, z(0)=0:
  z' ≤ (W √B - ν κ_H) z + U_L √B.
If α = W √B - ν κ_H > 0 and β = U_L √B,
  z(t) ≤ (β/α) (e^{α t} - 1),
  Ω(t) ≤ B + K_full² z(t)².

This improves the C-S-0002-SHORT horizon over L-0011.

FINITE Galerkin ONLY. Not continuum. Not Clay. Evidence: N7.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2
from ns_exploration.conjectures.l0009_twoscale import shell_mode_count
from ns_exploration.conjectures.l0010_duhamel_short import kappa_H_min
from ns_exploration.conjectures.l0011_sharpened_short import embedding_U, embedding_W
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def cancellation_z(
    t: float,
    U_L: float,
    W: float,
    B: float,
    nu: float,
    kappa_H: float,
) -> float:
    """Upper bound for √E_H from the linear comparison ODE."""
    if t <= 0:
        return 0.0
    beta = U_L * math.sqrt(B)
    alpha = W * math.sqrt(B) - nu * kappa_H
    if abs(alpha) < 1e-15:
        return beta * t  # z' ≤ β
    if alpha < 0:
        # z' ≤ -|α| z + β → z ≤ (β/|α|)(1 - e^{-|α|t}) ≤ β/|α|
        return (beta / (-alpha)) * (1.0 - math.exp(alpha * t))
    return (beta / alpha) * (math.exp(alpha * t) - 1.0)


def omega_cancellation(
    t: float,
    B: float,
    K_full: float,
    U_L: float,
    W: float,
    nu: float,
    kappa_H: float,
) -> float:
    z = cancellation_z(t, U_L, W, B, nu, kappa_H)
    return B + (K_full**2) * (z**2)


def max_T_cancellation(
    B: float,
    K_full: float,
    U_L: float,
    W: float,
    R: float,
    nu: float,
    kappa_H: float,
    T_hi: float = 0.05,
) -> float:
    if R <= B:
        return 0.0
    lo, hi = 0.0, T_hi
    if omega_cancellation(hi, B, K_full, U_L, W, nu, kappa_H) <= R:
        return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if omega_cancellation(mid, B, K_full, U_L, W, nu, kappa_H) <= R:
            lo = mid
        else:
            hi = mid
    return lo


@dataclass
class GalerkinBoundL0012:
    lemma_id: str = "L-0012"
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
    alpha: float = 0.0
    beta: float = 0.0
    cs0002_M: float = 48.15928260103266
    T_star: float = 0.0
    T_L0011_ref: float = 0.0012404305239949198
    improvement_vs_L0011: float = 0.0
    Omega_at_T_star: float = 0.0
    Omega_at_0_02: float = 0.0
    clay_implication: str = (
        "None. Finite Galerkin short-time bound via NS cancellations; "
        "does not prove C-S-0002 at T=0.02."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0012(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    nu: float = 0.1,
    cs0002_M: float = 48.15928260103266,
) -> GalerkinBoundL0012:
    if ic_kind != "linf":
        raise ValueError("L-0012 packaged for linf C-S class")

    K_ic2 = linf_shell_K2(k_ic)
    K2_full, M_full = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = shell_mode_count(n, k_ic, ic_kind)
    M_H = M_full - M_L
    kappa = float(kappa_H_min(n, k_ic, ic_kind) or 0.0)
    B = float(K_ic2) * E0
    U_L = embedding_U(M_L, E0, True)
    W = embedding_W(M_H, True)
    beta = U_L * math.sqrt(B)
    alpha = W * math.sqrt(B) - nu * kappa

    T_star = max_T_cancellation(B, K_full, U_L, W, cs0002_M, nu, kappa)
    om_star = omega_cancellation(T_star, B, K_full, U_L, W, nu, kappa)
    om_02 = omega_cancellation(0.02, B, K_full, U_L, W, nu, kappa)
    improv = T_star / 0.0012404305239949198

    notes = (
        f"NS cancellations: z' <= alpha z + beta with alpha={alpha:.6g}, beta={beta:.6g}. "
        f"T*={T_star:.8f} for R={cs0002_M:.4g} ({improv:.3g}x vs L-0011). "
        f"Omega(0.02)={om_02:.6g} (comparison ODE, may be huge)."
    )
    return GalerkinBoundL0012(
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
        alpha=alpha,
        beta=beta,
        cs0002_M=cs0002_M,
        T_star=T_star,
        improvement_vs_L0011=improv,
        Omega_at_T_star=om_star,
        Omega_at_0_02=om_02,
        notes=notes,
    )


def save_lemma_l0012(bound: GalerkinBoundL0012, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
