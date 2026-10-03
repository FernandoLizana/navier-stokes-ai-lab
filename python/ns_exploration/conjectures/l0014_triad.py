"""
L-0014: Shell→high triad multiplicity bound for (u_L·∇)u_L.

For û supported on the IC shell S_L, the convective coefficients satisfy
  |N̂_k| ≤ Σ_{p+q=k} |q| ||û_p|| ||û_q|| = (A ∗ B)_k,
with A_p = ||û_p||, B_q = |q| ||û_q||. Writing R_H = max_{k∈H} #{(p,q)∈S_L² : p+q=k},
  ||P_H N||_ℓ₂ ≤ √R_H ||A||_ℓ₂ ||B||_ℓ₂ = √R_H √(2E) √(2Ω) = 2 √(R_H E Ω).

Hence, under E≤E0 and Ω_L≤B,
  ||P_H (u_L·∇)u_L||_ℓ₂ ≤ C_△ √(E0 B),   C_△ := 2 √R_H,
and the L-0012 comparison ODE applies with the effective embedding
  U_eff = C_△ √(E0 / 2) = √(R_H) √(2 E0)
in place of the vector-CS U_L = √(M_L) √(2 E0). Since R_H ≤ M_L this is sharper.

The H×L remainder still uses W from L-0013.

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
    embedding_U_vector,
    embedding_W_vector,
    nonzero_high_mode_count,
    nonzero_shell_mode_count,
)
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def shell_high_triad_Rmax(n: int, k_ic: int, ic_kind: str = "linf") -> int:
    """Max number of shell+shell → high representations p+q=k for k in H."""
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
    ctr: Counter[tuple[int, int, int]] = Counter()
    for p in modes:
        for q in modes:
            s = (p[0] + q[0], p[1] + q[1], p[2] + q[2])
            if s in hset:
                ctr[s] += 1
    return int(max(ctr.values())) if ctr else 0


def triad_C_delta(R_H: int) -> float:
    """C_△ = 2 √R_H so that ||P_H N|| ≤ C_△ √(E Ω)."""
    return 2.0 * math.sqrt(float(R_H))


def effective_U_from_triads(R_H: int, E0: float) -> float:
    """U_eff = √R_H √(2 E0) for the L-0012 β = U_eff √B convention."""
    return math.sqrt(float(R_H)) * math.sqrt(2.0 * E0)


@dataclass
class GalerkinBoundL0014:
    lemma_id: str = "L-0014"
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
    R_H: int = 0
    C_delta: float = 0.0
    kappa_H: float = 0.0
    E0: float = 0.5
    nu: float = 0.1
    B: float = 0.0
    U_eff: float = 0.0
    U_L_L0013: float = 0.0
    W: float = 0.0
    alpha: float = 0.0
    beta: float = 0.0
    cs0002_M: float = 48.15928260103266
    T_star: float = 0.0
    T_L0013_ref: float = 0.0032350612481493943
    improvement_vs_L0013: float = 0.0
    Omega_at_T_star: float = 0.0
    Omega_at_0_02: float = 0.0
    clay_implication: str = (
        "None. Finite Galerkin short-time bound via shell→high triad count; "
        "does not prove C-S-0002 at T=0.02."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0014(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    nu: float = 0.1,
    cs0002_M: float = 48.15928260103266,
) -> GalerkinBoundL0014:
    if ic_kind != "linf":
        raise ValueError("L-0014 packaged for linf C-S class")

    K_ic2 = linf_shell_K2(k_ic)
    K2_full, _ = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = nonzero_shell_mode_count(n, k_ic, ic_kind)
    M_H = nonzero_high_mode_count(n, k_ic, ic_kind)
    R_H = shell_high_triad_Rmax(n, k_ic, ic_kind)
    C_delta = triad_C_delta(R_H)
    kappa = float(kappa_H_min(n, k_ic, ic_kind) or 0.0)
    B = float(K_ic2) * E0
    U_eff = effective_U_from_triads(R_H, E0)
    U_old = embedding_U_vector(M_L, E0)
    W = embedding_W_vector(M_H)

    beta = U_eff * math.sqrt(B)
    alpha = W * math.sqrt(B) - nu * kappa
    T_star = max_T_cancellation(B, K_full, U_eff, W, cs0002_M, nu, kappa)
    om_star = omega_cancellation(T_star, B, K_full, U_eff, W, nu, kappa)
    om_02 = omega_cancellation(0.02, B, K_full, U_eff, W, nu, kappa)
    improv = T_star / 0.0032350612481493943

    notes = (
        f"Triad R_H={R_H}, C_△={C_delta:.6g}, U_eff={U_eff:.6g} (L-0013 U_L={U_old:.6g}). "
        f"alpha={alpha:.6g}, beta={beta:.6g}. "
        f"T*={T_star:.8f} for R={cs0002_M:.4g} ({improv:.3g}x vs L-0013). "
        f"Omega(0.02)={om_02:.6g}."
    )
    return GalerkinBoundL0014(
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        R_H=R_H,
        C_delta=C_delta,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        B=B,
        U_eff=U_eff,
        U_L_L0013=U_old,
        W=W,
        alpha=alpha,
        beta=beta,
        cs0002_M=cs0002_M,
        T_star=T_star,
        improvement_vs_L0013=improv,
        Omega_at_T_star=om_star,
        Omega_at_0_02=om_02,
        notes=notes,
    )


def save_lemma_l0014(bound: GalerkinBoundL0014, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
