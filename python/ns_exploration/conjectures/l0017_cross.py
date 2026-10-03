"""
L-0017: Radial Young bound for the H×L cross term in the E_H ODE.

After NS cancellations (L-0012), the remainder involves
  ⟨u_H, (u_H·∇)u_L⟩.
Fourier: ((u_H·∇)u_L)_k = i Σ_{p∈H, q∈L, p+q=k} (û_p·q) û_q.
Young's inequality ‖A∗B‖₂ ≤ ‖A‖₂‖B‖₁ with A=‖û‖ on H and
B_q = |q|‖û_q‖ on L yields
  ‖(u_H·∇)u_L‖_ℓ₂ ≤ √(2 E_H) Σ_{q∈L} |q| ‖û_q‖.
The L-0016 radial estimate gives Σ_{q∈L} |q|‖û_q‖ ≤ √(2) √ρ_★ √E_L,
so with E_L≤E0
  ⟨u_H,(u_H·∇)u_L⟩ ≤ 2√2 √(ρ_★ E0) E_H.
In the L-0012 z=√E_H convention this replaces W√B by
  γ_cross = √(2 ρ_★ E0) = √ρ_★   (at E0=1/2),
and we take
  α = min(W √B, γ_cross) − ν κ_H
with U_eff / β from L-0016.

For N=16 ℓ∞≤4: γ_cross≈37.31 ≪ W√B≈170, T*≈0.01606.

FINITE Galerkin ONLY. Not continuum. Not Clay. Evidence: N7 (+ N5 cert).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2
from ns_exploration.conjectures.l0010_duhamel_short import kappa_H_min
from ns_exploration.conjectures.l0012_cancellation import omega_cancellation
from ns_exploration.conjectures.l0013_embedding import (
    embedding_W_vector,
    nonzero_high_mode_count,
    nonzero_shell_mode_count,
)
from ns_exploration.conjectures.l0016_radial import (
    Nmax_from_rho,
    effective_U_from_rho,
    shell_radial_rho_star,
)
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def gamma_cross_from_rho(rho_star: float, E0: float) -> float:
    """z-coefficient γ for the H×L term: z' ≤ … + γ z + …"""
    return math.sqrt(2.0 * float(rho_star) * E0)


def omega_cancellation_cross(
    t: float,
    B: float,
    K_full: float,
    U_eff: float,
    gamma_cross: float,
    W: float,
    nu: float,
    kappa_H: float,
) -> float:
    """Ω bound with α = min(W√B, γ_cross) − νκ_H, β = U_eff √B."""
    if t <= 0:
        return B
    beta = U_eff * math.sqrt(B)
    alpha = min(W * math.sqrt(B), gamma_cross) - nu * kappa_H
    if abs(alpha) < 1e-15:
        z = beta * t
    elif alpha < 0:
        z = (beta / (-alpha)) * (1.0 - math.exp(alpha * t))
    else:
        z = (beta / alpha) * (math.exp(alpha * t) - 1.0)
    return B + (K_full**2) * (z**2)


def max_T_cancellation_cross(
    B: float,
    K_full: float,
    U_eff: float,
    gamma_cross: float,
    W: float,
    R: float,
    nu: float,
    kappa_H: float,
    T_hi: float = 0.05,
) -> float:
    if R <= B:
        return 0.0
    lo, hi = 0.0, T_hi
    if omega_cancellation_cross(hi, B, K_full, U_eff, gamma_cross, W, nu, kappa_H) <= R:
        return hi
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if omega_cancellation_cross(mid, B, K_full, U_eff, gamma_cross, W, nu, kappa_H) <= R:
            lo = mid
        else:
            hi = mid
    return lo


@dataclass
class GalerkinBoundL0017:
    lemma_id: str = "L-0017"
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
    rho_star: float = 0.0
    r_star: int = 0
    m_star: int = 0
    N_max: float = 0.0
    kappa_H: float = 0.0
    E0: float = 0.5
    nu: float = 0.1
    B: float = 0.0
    U_eff: float = 0.0
    W: float = 0.0
    gamma_cross: float = 0.0
    W_sqrt_B: float = 0.0
    alpha: float = 0.0
    beta: float = 0.0
    cs0002_M: float = 48.15928260103266
    T_star: float = 0.0
    T_L0016_ref: float = 0.009115604958434453
    improvement_vs_L0016: float = 0.0
    Omega_at_T_star: float = 0.0
    Omega_at_0_02: float = 0.0
    clay_implication: str = (
        "None. Finite Galerkin short-time bound via radial H×L cross; "
        "does not prove C-S-0002 at T=0.02."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0017(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    nu: float = 0.1,
    cs0002_M: float = 48.15928260103266,
) -> GalerkinBoundL0017:
    if ic_kind != "linf":
        raise ValueError("L-0017 packaged for linf C-S class")

    K_ic2 = linf_shell_K2(k_ic)
    K2_full, _ = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = nonzero_shell_mode_count(n, k_ic, ic_kind)
    M_H = nonzero_high_mode_count(n, k_ic, ic_kind)
    rho, r_star, m_star = shell_radial_rho_star(n, k_ic, ic_kind)
    N_max = Nmax_from_rho(rho, E0)
    kappa = float(kappa_H_min(n, k_ic, ic_kind) or 0.0)
    B = float(K_ic2) * E0
    U_eff = effective_U_from_rho(rho, E0, K_ic2)
    W = embedding_W_vector(M_H)
    g_cross = gamma_cross_from_rho(rho, E0)
    W_sqrt_B = W * math.sqrt(B)

    beta = U_eff * math.sqrt(B)
    alpha = min(W_sqrt_B, g_cross) - nu * kappa
    T_star = max_T_cancellation_cross(
        B, K_full, U_eff, g_cross, W, cs0002_M, nu, kappa
    )
    om_star = omega_cancellation_cross(
        T_star, B, K_full, U_eff, g_cross, W, nu, kappa
    )
    om_02 = omega_cancellation_cross(
        0.02, B, K_full, U_eff, g_cross, W, nu, kappa
    )
    improv = T_star / 0.009115604958434453

    notes = (
        f"γ_cross={g_cross:.6g} vs W√B={W_sqrt_B:.6g}; U_eff={U_eff:.6g} (L-0016). "
        f"alpha={alpha:.6g}, beta={beta:.6g}. "
        f"T*={T_star:.8f} for R={cs0002_M:.4g} ({improv:.3g}x vs L-0016). "
        f"Omega(0.02)={om_02:.6g}."
    )
    return GalerkinBoundL0017(
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        rho_star=rho,
        r_star=r_star,
        m_star=m_star,
        N_max=N_max,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        B=B,
        U_eff=U_eff,
        W=W,
        gamma_cross=g_cross,
        W_sqrt_B=W_sqrt_B,
        alpha=alpha,
        beta=beta,
        cs0002_M=cs0002_M,
        T_star=T_star,
        improvement_vs_L0016=improv,
        Omega_at_T_star=om_star,
        Omega_at_0_02=om_02,
        notes=notes,
    )


def save_lemma_l0017(bound: GalerkinBoundL0017, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
