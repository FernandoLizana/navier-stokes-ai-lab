"""
L-0015: Energy-only shell→high bound via weighted triad row sums.

From Cauchy–Schwarz on each high mode,
  |N̂_k|² ≤ (Σ_{p+q=k} ||û_p||² |q|²) (Σ_{p+q=k} ||û_q||²)
         ≤ (Σ_{p+q=k} ||û_p||² |q|²) · (2E).
Summing k∈H and rearranging,
  ||P_H N||_ℓ₂² ≤ 2E · Σ_p ||û_p||² (Σ_{q: p+q∈H} |q|²)
                 ≤ 4 E² S_max,
where
  S_max := max_{p∈S_L} Σ_{q∈S_L : p+q∈H} |q|².
Hence under E≤E0,
  ||P_H (u_L·∇)u_L||_ℓ₂ ≤ 2 E0 √S_max =: N_max.

In the L-0012 ODE convention (dE production ≤ 2 U √B √E_H),
  U_eff = N_max / (√(2) √B) = √S_max · √(2 E0 / K_IC²).

For N=16, |k|_∞≤4: S_max=6132, U_eff≈11.30 (vs L-0014 U_eff=18).

FINITE Galerkin ONLY. Not continuum. Not Clay. Evidence: N7 (+ N5 cert).
"""

from __future__ import annotations

import json
import math
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
from ns_exploration.conjectures.l0014_triad import effective_U_from_triads, shell_high_triad_Rmax
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def shell_high_Smax(n: int, k_ic: int, ic_kind: str = "linf") -> float:
    """S_max = max_p Σ_{q: p+q∈H} |q|² over shell modes p."""
    if ic_kind != "linf":
        raise ValueError(ic_kind)
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    linf = np.maximum(np.maximum(np.abs(kx), np.abs(ky)), np.abs(kz))
    zero = (kx == 0) & (ky == 0) & (kz == 0)
    shell = mask & (linf <= float(k_ic) + 1e-12) & ~zero
    high = mask & (linf > float(k_ic) + 1e-12)
    modes = list(
        zip(
            kx[shell].astype(int),
            ky[shell].astype(int),
            kz[shell].astype(int),
        )
    )
    hset = set(
        zip(
            kx[high].astype(int),
            ky[high].astype(int),
            kz[high].astype(int),
        )
    )
    s_max = 0.0
    for p in modes:
        s = 0.0
        for q in modes:
            if (p[0] + q[0], p[1] + q[1], p[2] + q[2]) in hset:
                s += float(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])
        if s > s_max:
            s_max = s
    return float(s_max)


def Nmax_from_Smax(S_max: float, E0: float) -> float:
    """||P_H N|| ≤ 2 E0 √S_max."""
    return 2.0 * E0 * math.sqrt(float(S_max))


def effective_U_from_Smax(S_max: float, E0: float, K_ic_squared: int) -> float:
    """U_eff so that β = U_eff √B matches N_max in the L-0012 convention."""
    B = float(K_ic_squared) * E0
    N_max = Nmax_from_Smax(S_max, E0)
    return N_max / (math.sqrt(2.0) * math.sqrt(B))


@dataclass
class GalerkinBoundL0015:
    lemma_id: str = "L-0015"
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
    S_max: float = 0.0
    N_max: float = 0.0
    R_H_L0014: int = 0
    kappa_H: float = 0.0
    E0: float = 0.5
    nu: float = 0.1
    B: float = 0.0
    U_eff: float = 0.0
    U_L0014: float = 0.0
    W: float = 0.0
    alpha: float = 0.0
    beta: float = 0.0
    cs0002_M: float = 48.15928260103266
    T_star: float = 0.0
    T_L0014_ref: float = 0.0043669005451596726
    improvement_vs_L0014: float = 0.0
    Omega_at_T_star: float = 0.0
    Omega_at_0_02: float = 0.0
    clay_implication: str = (
        "None. Finite Galerkin short-time bound via S_max triad weights; "
        "does not prove C-S-0002 at T=0.02."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0015(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    nu: float = 0.1,
    cs0002_M: float = 48.15928260103266,
) -> GalerkinBoundL0015:
    if ic_kind != "linf":
        raise ValueError("L-0015 packaged for linf C-S class")

    K_ic2 = linf_shell_K2(k_ic)
    K2_full, _ = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = nonzero_shell_mode_count(n, k_ic, ic_kind)
    M_H = nonzero_high_mode_count(n, k_ic, ic_kind)
    S_max = shell_high_Smax(n, k_ic, ic_kind)
    R_H = shell_high_triad_Rmax(n, k_ic, ic_kind)
    N_max = Nmax_from_Smax(S_max, E0)
    kappa = float(kappa_H_min(n, k_ic, ic_kind) or 0.0)
    B = float(K_ic2) * E0
    U_eff = effective_U_from_Smax(S_max, E0, K_ic2)
    U_old = effective_U_from_triads(R_H, E0)
    W = embedding_W_vector(M_H)

    beta = U_eff * math.sqrt(B)
    alpha = W * math.sqrt(B) - nu * kappa
    T_star = max_T_cancellation(B, K_full, U_eff, W, cs0002_M, nu, kappa)
    om_star = omega_cancellation(T_star, B, K_full, U_eff, W, nu, kappa)
    om_02 = omega_cancellation(0.02, B, K_full, U_eff, W, nu, kappa)
    improv = T_star / 0.0043669005451596726

    notes = (
        f"S_max={S_max:.6g}, N_max={N_max:.6g}, U_eff={U_eff:.6g} "
        f"(L-0014 U_eff={U_old:.6g} from R_H={R_H}). "
        f"alpha={alpha:.6g}, beta={beta:.6g}. "
        f"T*={T_star:.8f} for R={cs0002_M:.4g} ({improv:.3g}x vs L-0014). "
        f"Omega(0.02)={om_02:.6g}."
    )
    return GalerkinBoundL0015(
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        S_max=S_max,
        N_max=N_max,
        R_H_L0014=R_H,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        B=B,
        U_eff=U_eff,
        U_L0014=U_old,
        W=W,
        alpha=alpha,
        beta=beta,
        cs0002_M=cs0002_M,
        T_star=T_star,
        improvement_vs_L0014=improv,
        Omega_at_T_star=om_star,
        Omega_at_0_02=om_02,
        notes=notes,
    )


def save_lemma_l0015(bound: GalerkinBoundL0015, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
