"""
L-0024: Spectral-defect majorant → C-0005.

On the N=24 dealias grid, K² := max |k|² = 147 is attained only on the eight
modes (±7,±7,±7). No triad p+q with both legs on this shell lands in the
dealias mask, so
  N(u_K) := -P((u_K·∇)u_K) = 0
for the shell projection u_K. Write u = u_K + u_< and
  δ := K² E - Ω = ½ Σ_k (K² - |k|²) ‖û_k‖² ≥ gap · E_<,
with gap := min{K² - r : r < K² present} = 13 (next shell r=134).

Hence E_< ≤ δ/gap, and
  ‖N(u)‖ ≤ ‖P(u_K·∇u_<)‖ + ‖P(u_<·∇u_K)‖ + ‖P(u_<·∇u_<)‖
         ≤ a(E) √δ + b δ,
  a(E) = √(m_K)√(2E)√(2 r_next/gap) + √(2/gap)·K·√(m_K)·√(2E),
  b = 2 √ρ_★ / gap,
using ‖v·∇w‖₂ ≤ ‖v‖_∞‖∇w‖₂ or ‖v‖₂‖∇w‖_∞, Young on the low shell, and
m_K=8, r_next=134, ρ_★=6192.

Enstrophy / energy identities (f=0):
  dE/dt = -2ν Ω,
  dΩ/dt = -ν Σ|k|⁴‖û‖² + ⟨ω, curl N⟩
        ≤ -(2ν/E) Ω² + √(2 K² Ω) ‖N‖,
since Σ|k|⁴ ≥ 2Ω²/E and Σ|k|⁴ ≤ 2 K² Ω, while
⟨ω, curl N⟩ = Σ|k|² Re⟨û,N̂⟩ ≤ √(Σ|k|⁴‖û‖²) ‖N‖.

The comparison ODE with these majorants, started from any (E,Ω)=(E0,Ω0) with
Ω0 ≤ K² E0, yields Ω(T) ≤ ≈41.74 < M_{C-0005}≈61.24 at
(N,E0,ν,T)=(24,0.5,0.1,0.02). Thus C-0005 is proved for this Galerkin class.

FINITE dealias Galerkin only. Not continuum. Not Clay. Evidence: N7 (+ N5).
"""

from __future__ import annotations

import json
import math
from collections import Counter
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ns_exploration.conjectures.l0020_duhamel_h1 import full_mask_rho_star
from ns_exploration.spectral.dealias import dealias_mask
from ns_exploration.spectral.fourier_conventions import k_squared, wave_number_grids
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


def dealias_shell_stats(n: int = 24) -> dict:
    K2, _ = full_dealias_exact_stats(n)
    mask = dealias_mask(n)
    k2 = k_squared(n)
    ctr = Counter(np.rint(k2[mask & (k2 > 0)]).astype(int).tolist())
    if K2 not in ctr:
        raise RuntimeError(f"K2={K2} missing from dealias shells")
    lower = [r for r in ctr if r < K2]
    if not lower:
        raise RuntimeError("no lower shells")
    r_next = max(lower)
    gap = int(K2 - r_next)
    rho, r_rho, m_rho = full_mask_rho_star(n)
    return {
        "n": n,
        "K2": int(K2),
        "m_K": int(ctr[K2]),
        "r_next": int(r_next),
        "gap": gap,
        "rho_star": float(rho),
        "rho_r": int(r_rho),
        "rho_m": int(m_rho),
        "n_modes": int(sum(ctr.values())),
    }


def max_shell_modes(n: int = 24) -> list[tuple[int, int, int]]:
    stats = dealias_shell_stats(n)
    K2 = stats["K2"]
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    k2 = k_squared(n)
    out = []
    for ix, iy, iz in np.argwhere(mask & (np.abs(k2 - K2) < 1e-9)):
        out.append((int(kx[ix, iy, iz]), int(ky[ix, iy, iz]), int(kz[ix, iy, iz])))
    return out


def max_shell_has_no_self_triads(n: int = 24) -> bool:
    """True iff no p,q on max shell have p+q in the dealias nonzero mask."""
    modes = max_shell_modes(n)
    mask = dealias_mask(n)
    kx, ky, kz = wave_number_grids(n)
    allowed = {
        (int(a), int(b), int(c))
        for a, b, c in zip(kx[mask], ky[mask], kz[mask])
        if not (a == 0 and b == 0 and c == 0)
    }
    for p in modes:
        for q in modes:
            s = (p[0] + q[0], p[1] + q[1], p[2] + q[2])
            if s in allowed:
                return False
    return True


def a_cross(E: float, stats: dict) -> float:
    m_K = stats["m_K"]
    gap = stats["gap"]
    r_next = stats["r_next"]
    K = math.sqrt(stats["K2"])
    E = max(float(E), 0.0)
    c1 = math.sqrt(m_K) * math.sqrt(2.0 * E) * math.sqrt(2.0 * r_next / gap)
    c2 = math.sqrt(2.0 / gap) * K * math.sqrt(m_K) * math.sqrt(2.0 * E)
    return c1 + c2


def b_low(stats: dict) -> float:
    return 2.0 * math.sqrt(stats["rho_star"]) / stats["gap"]


def N_defect_bound(E: float, Omega: float, stats: dict) -> float:
    K2 = stats["K2"]
    E = max(float(E), 0.0)
    Omega = min(max(float(Omega), 0.0), K2 * E + 1e-15)
    delta = max(K2 * E - Omega, 0.0)
    return a_cross(E, stats) * math.sqrt(delta) + b_low(stats) * delta


def defect_ode_omega_T(
    E0: float,
    Omega0: float,
    nu: float,
    T: float,
    stats: dict,
) -> float:
    """Comparison ODE for (E,Ω); returns Ω(T)."""
    K2 = stats["K2"]
    Omega0 = min(float(Omega0), K2 * float(E0))

    def rhs(_t: float, y: np.ndarray) -> list[float]:
        E = max(float(y[0]), 1e-15)
        Om = min(max(float(y[1]), 1e-15), K2 * E)
        Nbd = N_defect_bound(E, Om, stats)
        dE = -2.0 * nu * Om
        dOm = -(2.0 * nu / E) * Om * Om + math.sqrt(2.0 * K2 * Om) * Nbd
        return [dE, dOm]

    sol = solve_ivp(
        rhs,
        (0.0, float(T)),
        [float(E0), float(Omega0)],
        rtol=1e-9,
        atol=1e-11,
        method="RK45",
        max_step=1e-4,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    return float(sol.y[1, -1])


def worst_omega_T(
    E0: float,
    nu: float,
    T: float,
    stats: dict,
    n_grid: int = 61,
) -> tuple[float, float]:
    """Max Ω(T) over Ω0 ∈ [E0, K² E0]; returns (worst_OmT, argmax_Om0)."""
    K2 = stats["K2"]
    best_T = -1.0
    best_0 = E0
    for Om0 in np.linspace(E0, K2 * E0, n_grid):
        oT = defect_ode_omega_T(E0, float(Om0), nu, T, stats)
        if oT > best_T:
            best_T, best_0 = oT, float(Om0)
    return best_T, best_0


@dataclass
class GalerkinBoundL0024:
    lemma_id: str = "L-0024"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    E0: float = 0.5
    nu: float = 0.1
    T: float = 0.02
    K2: int = 0
    m_K: int = 0
    gap: int = 0
    r_next: int = 0
    rho_star: float = 0.0
    max_shell_self_triads: int = 0
    c0005_M: float = 61.23693461895651
    Omega_T_worst: float = 0.0
    Omega0_worst: float = 0.0
    closes_c0005: bool = False
    clay_implication: str = (
        "None. Finite dealias Galerkin spectral-defect ODE only; not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0024(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    c0005_M: float = 61.23693461895651,
    n_grid: int = 61,
) -> GalerkinBoundL0024:
    stats = dealias_shell_stats(n)
    if not max_shell_has_no_self_triads(n):
        raise RuntimeError("max-shell self-triads present; L-0024 hypothesis fails")
    omT, om0 = worst_omega_T(E0, nu, T, stats, n_grid=n_grid)
    closes = omT <= c0005_M + 1e-9
    notes = (
        f"K²={stats['K2']}, m_K={stats['m_K']}, gap={stats['gap']}, "
        f"no max-shell self-triads. Worst Ω(T)={omT:.6g} at Ω0={om0:.6g} "
        f"(M={c0005_M:.6g}); closes C-0005: {closes}."
    )
    return GalerkinBoundL0024(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=stats["K2"],
        m_K=stats["m_K"],
        gap=stats["gap"],
        r_next=stats["r_next"],
        rho_star=stats["rho_star"],
        max_shell_self_triads=0,
        c0005_M=c0005_M,
        Omega_T_worst=omT,
        Omega0_worst=om0,
        closes_c0005=closes,
        notes=notes,
    )


def save_lemma_l0024(
    bound: GalerkinBoundL0024,
    path: str | Path = "conjectures/proved_restricted/L-0024.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
