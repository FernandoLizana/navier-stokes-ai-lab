"""
L-0013: Corrected vector Fourier embedding + L-0012 cancellation ODE.

Pointwise Cauchy–Schwarz on Fourier modes (NS-MRL conventions):
  |u(x)| = |Σ_k û_k e^{ik·x}| ≤ Σ_k ||û_k||₂ ≤ √M √(Σ_k ||û_k||₂²) = √M √(2E),
where the sum is over nonzero lattice modes in the support (k=0 carries no
mean-zero velocity). Hence
  U_L = √(M_L) √(2 E0),   ||u_H||_∞ ≤ W √E_H with W = √(2 M_H).

This replaces the looser L-0011 factors √(2 M)√(2E) and √(4 M_H), which
inserted an extraneous √2 relative to the vector CS bound. Divergence-free
constraints only shrink the feasible set; they do not enlarge this upper bound.

The short-time ODE is exactly L-0012 (NS cancellations for E_H).

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
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def nonzero_shell_mode_count(n: int, k_ic: int, ic_kind: str = "linf") -> int:
    """Dealiased IC-shell modes excluding k=0."""
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    zero = (kx == 0) & (ky == 0) & (kz == 0)
    if ic_kind == "linf":
        linf = np.maximum(np.maximum(np.abs(kx), np.abs(ky)), np.abs(kz))
        sm = mask & (linf <= float(k_ic) + 1e-12) & ~zero
    else:
        raise ValueError(ic_kind)
    return int(sm.sum())


def nonzero_high_mode_count(n: int, k_ic: int, ic_kind: str = "linf") -> int:
    """Dealiased complement of the IC shell (automatically excludes k=0)."""
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    if ic_kind == "linf":
        linf = np.maximum(np.maximum(np.abs(kx), np.abs(ky)), np.abs(kz))
        low = mask & (linf <= float(k_ic) + 1e-12)
    else:
        raise ValueError(ic_kind)
    high = mask & ~low
    return int(high.sum())


def embedding_U_vector(M: int, E0: float) -> float:
    """||u||_∞ ≤ √M √(2E) for energy ≤ E0 on M nonzero modes."""
    if M <= 0:
        return 0.0
    return math.sqrt(float(M)) * math.sqrt(2.0 * E0)


def embedding_W_vector(M_H: int) -> float:
    """||u_H||_∞ ≤ √(2 M_H) √E_H."""
    if M_H <= 0:
        return 0.0
    return math.sqrt(2.0 * float(M_H))


@dataclass
class GalerkinBoundL0013:
    lemma_id: str = "L-0013"
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
    U_L_L0012: float = 0.0
    W_L0012: float = 0.0
    alpha: float = 0.0
    beta: float = 0.0
    cs0002_M: float = 48.15928260103266
    T_star: float = 0.0
    T_L0012_ref: float = 0.0022840472345178947
    improvement_vs_L0012: float = 0.0
    Omega_at_T_star: float = 0.0
    Omega_at_0_02: float = 0.0
    clay_implication: str = (
        "None. Finite Galerkin short-time bound; corrected embedding + L-0012 ODE. "
        "Does not prove C-S-0002 at T=0.02."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0013(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    nu: float = 0.1,
    cs0002_M: float = 48.15928260103266,
) -> GalerkinBoundL0013:
    if ic_kind != "linf":
        raise ValueError("L-0013 packaged for linf C-S class")

    K_ic2 = linf_shell_K2(k_ic)
    K2_full, _M_full = full_dealias_exact_stats(n)
    K_full = math.sqrt(float(K2_full))
    M_L = nonzero_shell_mode_count(n, k_ic, ic_kind)
    M_H = nonzero_high_mode_count(n, k_ic, ic_kind)
    kappa = float(kappa_H_min(n, k_ic, ic_kind) or 0.0)
    B = float(K_ic2) * E0
    U_L = embedding_U_vector(M_L, E0)
    W = embedding_W_vector(M_H)
    # L-0012 reference factors (including zero in shell count as historically used)
    from ns_exploration.conjectures.l0009_twoscale import shell_mode_count
    from ns_exploration.conjectures.l0011_sharpened_short import embedding_U, embedding_W

    M_L_old = shell_mode_count(n, k_ic, ic_kind)
    M_H_old = _M_full - M_L_old
    U_old = embedding_U(M_L_old, E0, True)
    W_old = embedding_W(M_H_old, True)

    beta = U_L * math.sqrt(B)
    alpha = W * math.sqrt(B) - nu * kappa
    T_star = max_T_cancellation(B, K_full, U_L, W, cs0002_M, nu, kappa)
    om_star = omega_cancellation(T_star, B, K_full, U_L, W, nu, kappa)
    om_02 = omega_cancellation(0.02, B, K_full, U_L, W, nu, kappa)
    improv = T_star / 0.0022840472345178947

    notes = (
        f"Vector CS embeddings: U_L={U_L:.6g} (was {U_old:.6g}), W={W:.6g} (was {W_old:.6g}); "
        f"M_L={M_L} excl. k=0 (was {M_L_old}). alpha={alpha:.6g}, beta={beta:.6g}. "
        f"T*={T_star:.8f} for R={cs0002_M:.4g} ({improv:.3g}x vs L-0012). "
        f"Omega(0.02)={om_02:.6g}."
    )
    return GalerkinBoundL0013(
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
        U_L_L0012=U_old,
        W_L0012=W_old,
        alpha=alpha,
        beta=beta,
        cs0002_M=cs0002_M,
        T_star=T_star,
        improvement_vs_L0012=improv,
        Omega_at_T_star=om_star,
        Omega_at_0_02=om_02,
        notes=notes,
    )


def save_lemma_l0013(bound: GalerkinBoundL0013, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
