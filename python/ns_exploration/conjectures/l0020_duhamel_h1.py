"""
L-0020: Stokes–Duhamel H¹ majorant for Ω(T) (toward C-0005).

Write the mild form on the dealias Galerkin space
  u(t) = e^{t ν Δ} u0 + ∫_0^t e^{(t-s) ν Δ} N(s) ds,
with N = P((u·∇)u). In Ḣ¹,
  ‖u(t)‖_Ḣ¹ ≤ ‖e^{t ν Δ} u0‖_Ḣ¹ + ∫_0^t ‖e^{(t-s) ν Δ} N(s)‖_Ḣ¹ ds.

Heat-kernel bounds on the dealias mask:
  ‖e^{t ν Δ} u0‖_Ḣ¹ ≤ √(2 S(t) E0),
  S(t) := max_r r e^{-2 ν r t},
  ‖e^{τ ν Δ} N‖_Ḣ¹ ≤ σ(τ) ‖N‖_ℓ₂,
  σ(τ) := max_r √r e^{-ν r τ}.

Hence if ‖N(s)‖_ℓ₂ ≤ N_* on [0,T],
  ‖u(T)‖_Ḣ¹ ≤ √(2 S(T) E0) + N_* ∫_0^T σ(τ) dτ,
  Ω(T) = (1/2) ‖u‖_Ḣ¹² ≤ (1/2)(√(2 S E0) + N_* I_σ)².

Unconditional: N_* = 2 E0 √ρ_★ (full-mask radial Young, L-0016 style) — too large
to beat the L-0018 envelope at these parameters.

Conditional for C-0005: need N_* ≤ N_crit where the bound equals M_C0005
(≈9.666 at N=24, T=0.02, E0=0.5, ν=0.1). Adversarial fields reach ‖N‖≳35
(see L-0022), so a uniform N_* instantiation cannot close all-IC C-0005;
unstructured random ICs alone can look like ‖N‖≲2 and are misleading.

FINITE Galerkin / dealias only. Not continuum. Not Clay. Evidence: N7 (+ N5).
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def dealias_nonzero_radii(n: int) -> list[int]:
    mask = dealias_mask(n)
    k2 = k_squared(n)
    rs = np.rint(k2[mask & (k2 > 0)]).astype(np.int64)
    return sorted({int(r) for r in rs.tolist() if int(r) > 0})


def full_mask_rho_star(n: int) -> tuple[float, int, int]:
    """ρ_★ = max_r (r m_r) on nonzero dealias modes."""
    mask = dealias_mask(n)
    k2 = k_squared(n)
    rs = np.rint(k2[mask & (k2 > 0)]).astype(int).tolist()
    ctr = Counter(rs)
    r_star, m_star, rho = 0, 0, 0.0
    for r, m in ctr.items():
        val = float(r * m)
        if val > rho:
            rho, r_star, m_star = val, int(r), int(m)
    return rho, r_star, m_star


def S_of_t(radii: list[int], nu: float, t: float) -> float:
    return max(float(r) * math.exp(-2.0 * nu * float(r) * t) for r in radii)


def sigma_of_tau(radii: list[int], nu: float, tau: float) -> float:
    return max(math.sqrt(float(r)) * math.exp(-nu * float(r) * tau) for r in radii)


def integrate_sigma(radii: list[int], nu: float, T: float, n_grid: int = 40001) -> float:
    if T <= 0:
        return 0.0
    u = np.linspace(0.0, T, n_grid)
    sig = np.array([sigma_of_tau(radii, nu, float(ui)) for ui in u], dtype=np.float64)
    return float(np.trapezoid(sig, u))


def young_Nstar(rho_star: float, E0: float) -> float:
    return 2.0 * E0 * math.sqrt(float(rho_star))


def omega_duhamel_H1(S_T: float, E0: float, N_star: float, I_sigma: float) -> float:
    h1 = math.sqrt(2.0 * S_T * E0) + float(N_star) * float(I_sigma)
    return 0.5 * h1 * h1


def Ncrit_for_M(S_T: float, E0: float, I_sigma: float, M: float) -> float:
    """Smallest N_* with omega_duhamel_H1 <= M (requires sqrt(2 S E0) <= sqrt(2M))."""
    st = math.sqrt(2.0 * S_T * E0)
    target = math.sqrt(2.0 * M)
    if target <= st or I_sigma <= 0:
        return 0.0
    return (target - st) / I_sigma


@dataclass
class GalerkinBoundL0020:
    lemma_id: str = "L-0020"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    E0: float = 0.5
    nu: float = 0.1
    T: float = 0.02
    K2: int = 0
    rho_star: float = 0.0
    r_star: int = 0
    m_star: int = 0
    S_T: float = 0.0
    Stokes_floor: float = 0.0
    I_sigma: float = 0.0
    Nstar_young: float = 0.0
    Omega_young_H1: float = 0.0
    envelope_cap: float = 0.0
    c0005_M: float = 61.23693461895651
    Ncrit_c0005: float = 0.0
    young_closes_c0005: bool = False
    young_beats_envelope: bool = False
    clay_implication: str = (
        "None. Finite dealias Galerkin Stokes–Duhamel H¹ bound; "
        "Young instantiation does not close C-0005; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0020(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    c0005_M: float = 61.23693461895651,
) -> GalerkinBoundL0020:
    radii = dealias_nonzero_radii(n)
    K2, _ = full_dealias_exact_stats(n)
    rho, r_star, m_star = full_mask_rho_star(n)
    S_T = S_of_t(radii, nu, T)
    I_sig = integrate_sigma(radii, nu, T)
    N_y = young_Nstar(rho, E0)
    om_y = omega_duhamel_H1(S_T, E0, N_y, I_sig)
    env = float(K2) * E0
    ncrit = Ncrit_for_M(S_T, E0, I_sig, c0005_M)
    notes = (
        f"S(T)={S_T:.6g}, Stokes_floor={S_T*E0:.6g}, I_σ={I_sig:.6g}. "
        f"Young N_*={N_y:.6g} (ρ_★={rho:.0f} at r={r_star},m={m_star}) ⇒ "
        f"Ω_H1={om_y:.6g} (envelope={env}). "
        f"C-0005 N_crit={ncrit:.6g}; Young closes C-0005: {om_y <= c0005_M}."
    )
    return GalerkinBoundL0020(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=int(K2),
        rho_star=rho,
        r_star=r_star,
        m_star=m_star,
        S_T=S_T,
        Stokes_floor=S_T * E0,
        I_sigma=I_sig,
        Nstar_young=N_y,
        Omega_young_H1=om_y,
        envelope_cap=env,
        c0005_M=c0005_M,
        Ncrit_c0005=ncrit,
        young_closes_c0005=om_y <= c0005_M + 1e-12,
        young_beats_envelope=om_y < env,
        notes=notes,
    )


def save_lemma_l0020(
    bound: GalerkinBoundL0020,
    path: str | Path = "conjectures/proved_restricted/L-0020.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
