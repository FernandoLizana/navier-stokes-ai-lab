"""
L-0016: Radial shell Young bound for shell→high convection.

Write N̂_k = Σ_{p+q=k} (û_p·q) û_q. With A_p=‖û_p‖, B_q=|q|‖û_q‖,
  ‖A ∗ B‖_ℓ₂ ≤ ‖A‖_ℓ₂ ‖B‖_ℓ₁,
so
  ‖P_H N‖_ℓ₂ ≤ √(2E) Σ_q |q| ‖û_q‖.

Group shell modes by r = |q|² with multiplicity m_r:
  Σ_q |q| ‖û_q‖ = Σ_r √r Σ_{|q|²=r} ‖û_q‖
                 ≤ Σ_r √(r m_r) √(2 E_r).

Maximizing over E_r≥0 with Σ E_r≤E0 puts all energy on the radius
maximizing r m_r:
  Σ_q |q| ‖û_q‖ ≤ √(2) √(ρ_★) √E0,   ρ_★ := max_r (r m_r).
Hence
  ‖P_H N‖ ≤ 2 E0 √ρ_★ =: N_max,
  U_eff = N_max / (√2 √B) = √ρ_★ · √(2 E0 / K_IC²).

For N=16, |k|_∞≤4: ρ_★=1392 (r=29, m_r=48), N_max≈37.31, U_eff≈5.39
(vs L-0015 U_eff≈11.30 from S_max).

FINITE Galerkin ONLY. Not continuum. Not Clay. Evidence: N7 (+ N5 cert).
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2
from ns_exploration.conjectures.l0010_duhamel_short import kappa_H_min
from ns_exploration.conjectures.l0012_cancellation import (
    max_T_cancellation,
    omega_cancellation,
)
from ns_exploration.conjectures.l0013_embedding import (
    embedding_W_vector,
    nonzero_high_mode_count,
    nonzero_shell_mode_count,
)
from ns_exploration.conjectures.l0015_smax import effective_U_from_Smax, shell_high_Smax
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def shell_radial_rho_star(n: int, k_ic: int, ic_kind: str = "linf") -> tuple[float, int, int]:
    """Return (ρ_★, r_★, m_★) with ρ_★ = max_r (r · m_r) on the IC shell."""
    if ic_kind != "linf":
        raise ValueError(ic_kind)
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    linf = np.maximum(np.maximum(np.abs(kx), np.abs(ky)), np.abs(kz))
    zero = (kx == 0) & (ky == 0) & (kz == 0)
    shell = mask & (linf <= float(k_ic) + 1e-12) & ~zero
    modes = list(
        zip(
            kx[shell].astype(int),
            ky[shell].astype(int),
            kz[shell].astype(int),
        )
    )
    ctr = Counter(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] for q in modes)
    r_star = 0
    m_star = 0
    rho = 0.0
    for r, m in ctr.items():
        val = float(r * m)
        if val > rho:
            rho = val
            r_star = int(r)
            m_star = int(m)
    return rho, r_star, m_star


def Nmax_from_rho(rho_star: float, E0: float) -> float:
    """‖P_H N‖ ≤ 2 E0 √ρ_★."""
    return 2.0 * E0 * math.sqrt(float(rho_star))


def effective_U_from_rho(rho_star: float, E0: float, K_ic_squared: int) -> float:
    B = float(K_ic_squared) * E0
    N_max = Nmax_from_rho(rho_star, E0)
    return N_max / (math.sqrt(2.0) * math.sqrt(B))


@dataclass
class GalerkinBoundL0016:
    lemma_id: str = "L-0016"
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
    S_max_L0015: float = 0.0
    kappa_H: float = 0.0
    E0: float = 0.5
    nu: float = 0.1
    B: float = 0.0
    U_eff: float = 0.0
    U_L0015: float = 0.0
    W: float = 0.0
    alpha: float = 0.0
    beta: float = 0.0
    cs0002_M: float = 48.15928260103266
    T_star: float = 0.0
    T_L0015_ref: float = 0.005967225327829871
    improvement_vs_L0015: float = 0.0
    Omega_at_T_star: float = 0.0
    Omega_at_0_02: float = 0.0
    clay_implication: str = (
        "None. Finite Galerkin short-time bound via radial Young estimate; "
        "does not prove C-S-0002 at T=0.02."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0016(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    nu: float = 0.1,
    cs0002_M: float = 48.15928260103266,
) -> GalerkinBoundL0016:
    if ic_kind != "linf":
        raise ValueError("L-0016 packaged for linf C-S class")

    K_ic2 = linf_shell_K2(k_ic)
    K2_full, _ = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = nonzero_shell_mode_count(n, k_ic, ic_kind)
    M_H = nonzero_high_mode_count(n, k_ic, ic_kind)
    rho, r_star, m_star = shell_radial_rho_star(n, k_ic, ic_kind)
    S_max = shell_high_Smax(n, k_ic, ic_kind)
    N_max = Nmax_from_rho(rho, E0)
    kappa = float(kappa_H_min(n, k_ic, ic_kind) or 0.0)
    B = float(K_ic2) * E0
    U_eff = effective_U_from_rho(rho, E0, K_ic2)
    U_old = effective_U_from_Smax(S_max, E0, K_ic2)
    W = embedding_W_vector(M_H)

    beta = U_eff * math.sqrt(B)
    alpha = W * math.sqrt(B) - nu * kappa
    T_star = max_T_cancellation(B, K_full, U_eff, W, cs0002_M, nu, kappa)
    om_star = omega_cancellation(T_star, B, K_full, U_eff, W, nu, kappa)
    om_02 = omega_cancellation(0.02, B, K_full, U_eff, W, nu, kappa)
    improv = T_star / 0.005967225327829871

    notes = (
        f"ρ_★={rho:.6g} at r={r_star} (m_r={m_star}); N_max={N_max:.6g}; "
        f"U_eff={U_eff:.6g} (L-0015 U_eff={U_old:.6g}). "
        f"alpha={alpha:.6g}, beta={beta:.6g}. "
        f"T*={T_star:.8f} for R={cs0002_M:.4g} ({improv:.3g}x vs L-0015). "
        f"Omega(0.02)={om_02:.6g}."
    )
    return GalerkinBoundL0016(
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
        S_max_L0015=S_max,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        B=B,
        U_eff=U_eff,
        U_L0015=U_old,
        W=W,
        alpha=alpha,
        beta=beta,
        cs0002_M=cs0002_M,
        T_star=T_star,
        improvement_vs_L0015=improv,
        Omega_at_T_star=om_star,
        Omega_at_0_02=om_02,
        notes=notes,
    )


def save_lemma_l0016(bound: GalerkinBoundL0016, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
