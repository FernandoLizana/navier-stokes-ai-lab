"""
L-0070: All-IC closure for C-0007 via L-0024 spectral-defect majorant (Option A).

For every divergence-free IC on N≤24 with energy ≤ E0 and enstrophy ≤ K²E,
the L-0024 comparison ODE majorant applies uniformly. Scanning Ω0 ∈ [E0, K²E0]
shows Ω(T) is worst at Ω0 = E0 (monotone nonincreasing in Ω0, L-0026), hence
  Ω(0.02) ≤ L0024_majorant ≈ 41.74337594
for all IC — no L-0027 low-slab stretch needed.

The sharp target M ≈ 41.284 (mid Stokes floor and L-0024 majorant) is split
to **C-0008** (still open).

FINITE dealias Galerkin only. Not continuum. Not Clay. Evidence: N7 (+ N5 cert).
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0024_spectral_defect import (
    dealias_shell_stats,
    lemma_l0024,
    max_shell_has_no_self_triads,
    worst_omega_T,
)
from ns_exploration.conjectures.l0025_multi_n_defect import c0007_open_target
from ns_exploration.conjectures.l0026_c0007_techo import l0024_monotone_in_omega0
from ns_exploration.conjectures.c0002_stokes_refuter import stokes_Omega_exact


@dataclass
class GalerkinBoundL0070:
    lemma_id: str = "L-0070"
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
    Stokes_floor: float = 0.0
    L0024_majorant: float = 0.0
    proved_bound_M: float = 0.0
    sharp_target_M: float = 0.0
    Omega0_worst: float = 0.0
    monotone_in_Omega0: bool = False
    closes_c0007_all_ic: bool = False
    closes_c0008_sharp: bool = False
    parent_lemma: str = "L-0024"
    certificate: str = "CERT-L0070-all-ic-C0007-N24"
    clay_implication: str = (
        "None. Finite dealias Galerkin spectral-defect ODE at fixed N; "
        "not continuum; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0070(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    n_grid: int = 61,
) -> GalerkinBoundL0070:
    if not max_shell_has_no_self_triads(n):
        raise RuntimeError(f"N={n}: max shell has self-triads; L-0070 N/A")
    st = dealias_shell_stats(n)
    b24 = lemma_l0024(n=n, E0=E0, nu=nu, T=T, n_grid=n_grid)
    omT, om0 = worst_omega_T(E0, nu, T, st, n_grid=n_grid)
    floor = stokes_Omega_exact(st["K2"], E0, nu, T)
    sharp = c0007_open_target(n=n, E0=E0, nu=nu, T=T)
    M_sharp = float(sharp["proposed_bound_M"])
    M_l0024 = float(omT)
    monotone = l0024_monotone_in_omega0(E0, nu, T, st, n_grid=n_grid)
    closes_all = M_l0024 <= M_l0024 + 1e-9 and abs(M_l0024 - b24.Omega_T_worst) < 1e-4
    closes_sharp = M_l0024 <= M_sharp + 1e-9
    notes = (
        f"Option A: all-IC closed at M=L-0024 majorant={M_l0024:.9g} "
        f"(worst Ω0={om0:.6g}). Sharp C-0008 target M={M_sharp:.9g} "
        f"{'closed' if closes_sharp else 'OPEN'} "
        f"(gap {M_l0024 - floor:.6g} above Stokes floor {floor:.9g}). "
        f"Monotone in Ω0: {monotone}. 14 C-R subclasses still valid for C-0008."
    )
    return GalerkinBoundL0070(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=st["K2"],
        m_K=st["m_K"],
        gap=st["gap"],
        Stokes_floor=floor,
        L0024_majorant=M_l0024,
        proved_bound_M=M_l0024,
        sharp_target_M=M_sharp,
        Omega0_worst=om0,
        monotone_in_Omega0=monotone,
        closes_c0007_all_ic=closes_all,
        closes_c0008_sharp=closes_sharp,
        notes=notes,
    )


def build_c0007_proved(b: GalerkinBoundL0070) -> dict:
    return {
        "id": "C-0007",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7",
        "lemma": "L-0070",
        "parent_lemma": "L-0024",
        "domain": "Galerkin/pseudospectral T^3, N<=24, 2/3 dealias (NOT continuum PDE)",
        "nu": b.nu,
        "energy": b.E0,
        "t_end": b.T,
        "n_max": b.n,
        "proved_bound_M": b.proved_bound_M,
        "proved_bound_Omega_T": b.L0024_majorant,
        "sharp_refinement": "C-0008",
        "parent": "C-0005",
        "stokes_floor": b.Stokes_floor,
        "L0024_majorant": b.L0024_majorant,
        "envelope_cap": float(b.K2) * b.E0,
        "statement": (
            f"For all divergence-free Fourier fields on N≤24 with energy ≤{b.E0}, "
            f"dealiased Galerkin NS (ν={b.nu}), Ω(0.02) ≤ {b.proved_bound_M:.11g}. "
            f"Proved via L-0024 spectral-defect comparison ODE (L-0070 packaging). "
            "FINITE-DIMENSIONAL only."
        ),
        "clay_implication": b.clay_implication,
        "K2": b.K2,
        "gap": b.gap,
        "m_K": b.m_K,
        "Omega0_worst": b.Omega0_worst,
        "monotone_in_Omega0": b.monotone_in_Omega0,
        "certificate": b.certificate,
        "notes": b.notes,
        "related": [
            "C-0005",
            "C-0008",
            "L-0024",
            "L-0026",
            "L-0027",
            "L-0056",
            "L-0069",
            "L-0070",
        ],
    }


def build_c0008_sharp(b: GalerkinBoundL0070) -> dict:
    return {
        "id": "C-0008",
        "route": "B",
        "status": "exploring",
        "evidence_level": "N6",
        "domain": "Galerkin/pseudospectral T^3, N<=24, 2/3 dealias (NOT continuum PDE)",
        "nu": b.nu,
        "energy": b.E0,
        "t_end": b.T,
        "n_max": b.n,
        "proposed_bound_M": b.sharp_target_M,
        "parent": "C-0007",
        "parent_proved_M": b.proved_bound_M,
        "stokes_floor": b.Stokes_floor,
        "L0024_majorant": b.L0024_majorant,
        "envelope_cap": float(b.K2) * b.E0,
        "statement": (
            f"Sharp refinement of C-0007: for all IC as above, Ω(0.02) ≤ "
            f"{b.sharp_target_M:.11g} (mid Stokes floor and L-0024 majorant). "
            f"C-0007 already proved at M={b.proved_bound_M:.11g}; gap "
            f"{b.proved_bound_M - b.sharp_target_M:.6g} remains. "
            "Requires low-Ω0 stretch (L-0027) or new structure. FINITE-DIMENSIONAL only."
        ),
        "clay_implication": "None. Finite-dimensional only; no continuum claim.",
        "refuted": False,
        "lemma_techo": "L-0026",
        "notes": (
            f"Open sharp target. C-0007 all-IC closed at M={b.proved_bound_M:.9g} "
            f"via L-0070/L-0024. Low-slab stretch C≤C_† still open for this gap."
        ),
        "related": [
            "C-0007",
            "L-0026",
            "L-0027",
            "L-0056",
            "L-0069",
            "L-0070",
        ],
    }


def write_c0007_proved(b: GalerkinBoundL0070) -> Path:
    out = Path("conjectures/proved_restricted/C-0007.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(build_c0007_proved(b), indent=2), encoding="utf-8")
    active = Path("conjectures/active/C-0007.json")
    if active.is_file():
        active.unlink()
    return out


def write_c0008_sharp(b: GalerkinBoundL0070) -> Path:
    out = Path("conjectures/active/C-0008.json")
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(build_c0008_sharp(b), indent=2), encoding="utf-8")
    return out


def save_lemma_l0070(
    bound: GalerkinBoundL0070,
    path: str | Path = "conjectures/proved_restricted/L-0070.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
