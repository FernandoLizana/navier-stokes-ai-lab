"""
L-0025: Spectral-defect majorant at general dealias N (extends L-0024).

Same argument as L-0024 whenever the maximal shell has no self-triads under
2/3 dealias (verified for N∈{12,16,24,32,48}: always the eight corners
(±q,±q,±q)). Instantiates C-0006 on N≤32.

FINITE dealias Galerkin only. Not continuum. Not Clay. Evidence: N7 (+ N5).
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.c0002_stokes_refuter import stokes_Omega_exact
from ns_exploration.conjectures.l0024_spectral_defect import (
    dealias_shell_stats,
    lemma_l0024,
    max_shell_has_no_self_triads,
    worst_omega_T,
)
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class GalerkinBoundL0025:
    lemma_id: str = "L-0025"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 32
    E0: float = 0.5
    nu: float = 0.1
    T: float = 0.02
    K2: int = 0
    m_K: int = 0
    gap: int = 0
    r_next: int = 0
    rho_star: float = 0.0
    Stokes_floor: float = 0.0
    envelope_cap: float = 0.0
    c0006_M: float = 0.0
    Omega_T_worst: float = 0.0
    Omega0_worst: float = 0.0
    closes_c0006: bool = False
    clay_implication: str = (
        "None. Finite dealias Galerkin spectral-defect ODE at fixed N; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0025(
    n: int = 32,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    cushion: float = 1.5,
    n_grid: int = 61,
) -> GalerkinBoundL0025:
    if not max_shell_has_no_self_triads(n):
        raise RuntimeError(f"N={n}: max shell has self-triads; L-0025 N/A")
    # Reuse L-0024 engine at this N
    b24 = lemma_l0024(n=n, E0=E0, nu=nu, T=T, c0005_M=1e100, n_grid=n_grid)
    K2, _ = full_dealias_exact_stats(n)
    floor = stokes_Omega_exact(K2, E0, nu, T)
    env = float(K2) * E0
    M6 = float(cushion) * floor
    closes = b24.Omega_T_worst <= M6 + 1e-9
    notes = (
        f"L-0024@N={n}: K²={b24.K2}, m_K={b24.m_K}, gap={b24.gap}, "
        f"Ω(T)_worst={b24.Omega_T_worst:.6g}, Stokes_floor={floor:.6g}, "
        f"C-0006 M={M6:.6g}; closes: {closes}."
    )
    return GalerkinBoundL0025(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=b24.K2,
        m_K=b24.m_K,
        gap=b24.gap,
        r_next=b24.r_next,
        rho_star=b24.rho_star,
        Stokes_floor=floor,
        envelope_cap=env,
        c0006_M=M6,
        Omega_T_worst=b24.Omega_T_worst,
        Omega0_worst=b24.Omega0_worst,
        closes_c0006=closes,
        notes=notes,
    )


def c0007_open_target(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
) -> dict:
    """
    Open gap below the L-0024 majorant but above the Stokes floor on N=24.
    M = mid(floor, OmT_worst) — requires a sharper-than-L-0024 argument.
    """
    st = dealias_shell_stats(n)
    floor = stokes_Omega_exact(st["K2"], E0, nu, T)
    omT, _ = worst_omega_T(E0, nu, T, st, n_grid=61)
    M = 0.5 * (floor + omT)
    return {
        "n": n,
        "Stokes_floor": floor,
        "L0024_majorant": omT,
        "proposed_bound_M": M,
        "envelope_cap": float(st["K2"]) * E0,
        "open_gap": omT - floor,
    }


def save_lemma_l0025(
    bound: GalerkinBoundL0025,
    path: str | Path = "conjectures/proved_restricted/L-0025.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
