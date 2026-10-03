"""
L-0068: Prove restricted subclass C-R-0014 (one-pol 12 shells incl. 24).

Extends C-R-0013; discovered by L-0064 greedy-second growth to r≤25.

FINITE Galerkin only. Not all-IC. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0043_onepol_shellblock import sym_onepol_C_shor
from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import SUPPORT_GREEDY_SECOND
from ns_exploration.conjectures.l0067_onepol_greedy_r24_sos import (
    SUPPORT_CR0014,
    POL_BRANCH,
    lemma_l0067,
)

SUPPORT_CR0013 = tuple(SUPPORT_GREEDY_SECOND)


@dataclass
class GalerkinBoundL0068:
    lemma_id: str = "L-0068"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    support_cr0014: list[int] | None = None
    pol_branch: str = POL_BRANCH
    support_cr0014_D: int = 0
    C_shor_sym_cr0014: float = 0.0
    C_ub_sos_hi: float | None = None
    closes_cr0014: bool = False
    extends_cr0013: bool = False
    closes_c0007_all_ic: bool = False
    clay_implication: str = "None. One-pol C-R-0014; not all-IC; not Clay."
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0068(n: int = 24, *, load_sos: bool = True) -> GalerkinBoundL0068:
    b7 = lemma_l0027(n=n, empirical=False)
    b6 = lemma_l0026(n=n, n_grid=21)
    Cd = b7.C_dagger
    ref = sym_onepol_C_shor(n, SUPPORT_CR0014, branch=POL_BRANCH)
    C = float(ref["C_shor_sym"])
    assert set(SUPPORT_CR0013).issubset(set(SUPPORT_CR0014))
    assert C <= Cd + 1e-9

    sos_hi: float | None = None
    if load_sos:
        try:
            sos_hi = lemma_l0067(n=n).C_ub_hi
            assert sos_hi < Cd - 1e-9
        except FileNotFoundError:
            pass

    return GalerkinBoundL0068(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        Omega_star=b6.Omega_star,
        support_cr0014=list(ref["radii"]),
        pol_branch=POL_BRANCH,
        support_cr0014_D=int(ref["D"]),
        C_shor_sym_cr0014=C,
        C_ub_sos_hi=sos_hi,
        closes_cr0014=True,
        extends_cr0013=True,
        notes=(
            f"C-R-0014: {len(SUPPORT_CR0014)} shells D={ref['D']} C_Shor={C:.4g}≤C_† "
            f"(margin≈{Cd - C:.4f}). Extends C-R-0013 by shell 24."
            + (f" SOS N5≈{sos_hi:.4f}." if sos_hi else "")
        ),
    )


def cr0014_record(bound: GalerkinBoundL0068) -> dict:
    sos = (
        f" SOS N5 C_ub_hi≈{bound.C_ub_sos_hi:.6f} (L-0067)."
        if bound.C_ub_sos_hi is not None
        else ""
    )
    return {
        "id": "C-R-0014",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7+N5" if bound.C_ub_sos_hi else "N7",
        "lemma": "L-0068+L-0026+L-0027" + ("+L-0067" if bound.C_ub_sos_hi else ""),
        "domain": (
            f"Pseudospectral T^3, N<=24, 2/3 dealias, one-pol branch={bound.pol_branch}, "
            f"shells |k|^2 in {list(bound.support_cr0014 or [])}. NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": bound.n,
        "proposed_bound_M": bound.c0007_M,
        "shells": list(bound.support_cr0014 or []),
        "pol_branch": bound.pol_branch,
        "C_shor_sym": bound.C_shor_sym_cr0014,
        "C_ub_sos_hi": bound.C_ub_sos_hi,
        "C_dagger": bound.C_dagger,
        "D": bound.support_cr0014_D,
        "Omega_star": bound.Omega_star,
        "parent_open": "C-0007",
        "extends": "C-R-0013",
        "statement": (
            f"For divergence-free fields on N≤{bound.n} with one-pol support on "
            f"shells {list(bound.support_cr0014 or [])}, dealiased Galerkin NS "
            f"(ν=0.1, E≤0.5) satisfies Ω(0.02)≤{bound.c0007_M}. "
            f"Shor C≤{bound.C_shor_sym_cr0014}≤C_† ⇒ L-0027; high slab L-0026.{sos} FINITE only."
        ),
        "clay_implication": bound.clay_implication,
        "notes": "L-0064 greedy-second growth; adds shell 24 to C-R-0013 witness.",
    }


def save_lemma_l0068(
    bound: GalerkinBoundL0068,
    path: str | Path = "conjectures/proved_restricted/L-0068.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0014(
    bound: GalerkinBoundL0068,
    path: str | Path = "conjectures/proved_restricted/C-R-0014.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cr0014_record(bound), indent=2), encoding="utf-8")
