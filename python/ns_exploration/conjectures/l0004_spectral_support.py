"""
L-0004: Tighter Galerkin bounds under explicit spectral support hypotheses.

Hypothesis H(S): velocity Fourier coefficients vanish outside a retained set S
of spatial lattice points (after dealias), with
  K_S = max_{k in S} |k|,   M_S = |S|.

Then the L-0001/L-0003 estimates hold with (K,M) <- (K_S, M_S).

Additionally, for a *given* field u_hat, define effective stats from energy:
  - sort modes by per-mode energy e_k = sum_i |û_k^i|^2
  - M_eff = min count needed to capture `energy_fraction` of total energy
  - K_eff = max |k| among those modes

Field-conditional bound (diagnostic): apply L-0003 with (K_eff, M_eff).
This is NOT valid for all fields unless H(S) is assumed.

Evidence: N7 for shell-restricted; N2 for field-conditional diagnostics.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0001_galerkin_bound import explicit_enstrophy_bound
from ns_exploration.conjectures.l0003_uniform_galerkin import explicit_uniform_bound
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared, wave_number_grids


@dataclass
class SpectralSupportStats:
    K_max: float
    M_total: int
    K_eff: float
    M_eff: int
    energy_fraction: float
    k_shell_max: int | None = None

    def as_dict(self) -> dict:
        return asdict(self)


def spectral_support_stats(
    u_hat: np.ndarray,
    energy_fraction: float = 0.99,
    dealias: bool = True,
    k_shell_max: int | None = None,
) -> SpectralSupportStats:
    n = u_hat.shape[-1]
    kx, ky, kz = wave_number_grids(n)
    k2 = k_squared(n)
    kr = np.sqrt(k2)
    mask = dealias_mask(n) if dealias else np.ones_like(k2, dtype=bool)
    if k_shell_max is not None:
        mask &= kr <= float(k_shell_max) + 1e-12

    e_mode = np.sum(np.abs(u_hat[:, mask]) ** 2, axis=0)
    flat_kr = kr[mask]
    M_total = int(e_mode.size)
    if M_total == 0:
        return SpectralSupportStats(0.0, 0, 0.0, 0, energy_fraction, k_shell_max)

    order = np.argsort(e_mode.ravel())[::-1]
    e_sorted = e_mode.ravel()[order]
    total = float(np.sum(e_sorted))
    if total <= 0:
        return SpectralSupportStats(
            float(np.max(flat_kr)), M_total, 0.0, 1, energy_fraction, k_shell_max
        )
    cum = np.cumsum(e_sorted) / total
    idx = int(np.searchsorted(cum, energy_fraction)) + 1
    idx = min(idx, M_total)
    active = order[:idx]
    K_eff = float(np.max(flat_kr.ravel()[active]))
    K_max = float(np.max(flat_kr))
    return SpectralSupportStats(
        K_max=K_max,
        M_total=M_total,
        K_eff=K_eff,
        M_eff=idx,
        energy_fraction=energy_fraction,
        k_shell_max=k_shell_max,
    )


@dataclass
class GalerkinBoundL0004:
    lemma_id: str = "L-0004"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    hypothesis: str = ""
    K: float = 0.0
    M: int = 0
    E0: float = 0.5
    nu: float = 0.1
    t: float = 0.02
    L0001_short: float = 0.0
    L0003_uniform: float = 0.0
    best_cap: float = 0.0
    clay_implication: str = "None. Support-restricted finite Galerkin only."
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def shell_a_priori_bounds(
    n: int,
    k_shell_max: int,
    E0: float = 0.5,
    t: float = 0.02,
    nu: float = 0.1,
    dealias: bool = True,
) -> GalerkinBoundL0004:
    """A priori bounds for all fields supported in |k| <= k_shell_max (and dealias)."""
    kx, ky, kz = wave_number_grids(n)
    k2 = k_squared(n)
    kr = np.sqrt(k2)
    mask = dealias_mask(n) if dealias else np.ones_like(k2, dtype=bool)
    mask &= kr <= float(k_shell_max) + 1e-12
    k_ret = kr[mask]
    K = float(np.max(k_ret)) if k_ret.size else 0.0
    M = int(np.count_nonzero(mask))
    b1 = explicit_enstrophy_bound(E0, t, K, M)
    b3 = explicit_uniform_bound(E0, K, M, nu=nu)
    best = min(b1.Omega_t_cap, b3.Omega_uniform_cap)
    return GalerkinBoundL0004(
        hypothesis=f"|k|<= {k_shell_max} shell on N={n} grid (dealias={dealias})",
        K=K,
        M=M,
        E0=E0,
        nu=nu,
        t=t,
        L0001_short=b1.Omega_t_cap,
        L0003_uniform=b3.Omega_uniform_cap,
        best_cap=best,
        notes=(
            f"Shell-restricted: K={K:.6g}, M={M}. "
            f"L-0001@t={t}: {b1.Omega_t_cap:.6e}; L-0003 uniform: {b3.Omega_uniform_cap:.6e}."
        ),
    )


def field_conditional_bounds(
    u_hat: np.ndarray,
    E0: float = 0.5,
    t: float = 0.02,
    nu: float = 0.1,
    energy_fraction: float = 0.99,
) -> GalerkinBoundL0004:
    """Diagnostic bounds using effective (K_eff, M_eff) — NOT universal."""
    st = spectral_support_stats(u_hat, energy_fraction=energy_fraction)
    b1 = explicit_enstrophy_bound(E0, t, st.K_eff, st.M_eff)
    b3 = explicit_uniform_bound(E0, st.K_eff, st.M_eff, nu=nu)
    best = min(b1.Omega_t_cap, b3.Omega_uniform_cap)
    return GalerkinBoundL0004(
        hypothesis=f"field-conditional M_eff={st.M_eff} at {energy_fraction*100:.1f}% energy",
        K=st.K_eff,
        M=st.M_eff,
        E0=E0,
        nu=nu,
        t=t,
        L0001_short=b1.Omega_t_cap,
        L0003_uniform=b3.Omega_uniform_cap,
        best_cap=best,
        evidence_level="N2",
        status="diagnostic",
        notes=f"K_eff={st.K_eff:.6g}, M_eff={st.M_eff}, M_total={st.M_total}. Not a universal bound.",
    )


def save_l0004(b: GalerkinBoundL0004, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(b.as_dict(), indent=2), encoding="utf-8")
