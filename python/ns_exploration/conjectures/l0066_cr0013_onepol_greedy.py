"""
L-0066: Prove restricted subclass C-R-0013 (one-pol greedy-second, L-0064/65).

Hermitian one-pol fields with support on 11 shells (branch=second) satisfy
C-0007 via Shor C≤C_† (N7). SOS N5 interval via L-0065 strengthens C bound.

FINITE Galerkin only. Not all-IC. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0043_onepol_shellblock import SUPPORT_CR0008, sym_onepol_C_shor
from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import (
    SUPPORT_GREEDY_SECOND,
    POL_BRANCH,
    lemma_l0065,
)

SUPPORT_CR0013 = tuple(SUPPORT_GREEDY_SECOND)


@dataclass
class GalerkinBoundL0066:
    lemma_id: str = "L-0066"
    route: str = "B"
    status: str = "proved_restricted"
    evidence_level: str = "N7"
    n: int = 24
    c0007_M: float = 0.0
    C_dagger: float = 0.0
    Omega_star: float = 0.0
    support_cr0013: list[int] | None = None
    pol_branch: str = POL_BRANCH
    support_cr0013_D: int = 0
    C_shor_sym_cr0013: float = 0.0
    C_ub_sos_hi: float | None = None
    closes_cr0013: bool = False
    extends_cr0008: bool = False
    closes_c0007_all_ic: bool = False
    clay_implication: str = (
        "None. One-pol greedy-second subclass C-R-0013; not all-IC; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0066(n: int = 24, *, load_sos: bool = True) -> GalerkinBoundL0066:
    b7 = lemma_l0027(n=n, empirical=False)
    b6 = lemma_l0026(n=n, n_grid=21)
    Cd = b7.C_dagger
    ref = sym_onepol_C_shor(n, SUPPORT_CR0013, branch=POL_BRANCH)
    C = float(ref["C_shor_sym"])
    assert set(SUPPORT_CR0008).issubset(set(SUPPORT_CR0013))
    assert C <= Cd + 1e-9
    assert C >= sym_onepol_C_shor(n, SUPPORT_CR0008, branch=POL_BRANCH)["C_shor_sym"] - 1e-6

    sos_hi: float | None = None
    if load_sos:
        try:
            sos_hi = lemma_l0065(n=n).C_ub_hi
            assert sos_hi < Cd - 1e-9
        except FileNotFoundError:
            pass

    notes = (
        f"C-R-0013 one-pol {POL_BRANCH}: {len(SUPPORT_CR0013)} shells D={ref['D']} "
        f"C_Shor={C:.4g}≤C_†. Extends C-R-0008 ({len(SUPPORT_CR0008)} shells)."
    )
    if sos_hi is not None:
        notes += f" SOS N5 (L-0065) C_ub_hi≈{sos_hi:.4f}."
    notes += " All-IC: False."

    return GalerkinBoundL0066(
        n=n,
        c0007_M=b7.c0007_M,
        C_dagger=Cd,
        Omega_star=b6.Omega_star,
        support_cr0013=list(ref["radii"]),
        pol_branch=POL_BRANCH,
        support_cr0013_D=int(ref["D"]),
        C_shor_sym_cr0013=C,
        C_ub_sos_hi=sos_hi,
        closes_cr0013=True,
        extends_cr0008=True,
        closes_c0007_all_ic=False,
        notes=notes,
    )


def cr0013_record(bound: GalerkinBoundL0066) -> dict:
    sos_note = (
        f" SOS N5 interval C_ub_hi≈{bound.C_ub_sos_hi:.6f} (L-0065)."
        if bound.C_ub_sos_hi is not None
        else ""
    )
    return {
        "id": "C-R-0013",
        "route": "B",
        "status": "proved_restricted",
        "evidence_level": "N7+N5" if bound.C_ub_sos_hi is not None else "N7",
        "lemma": "L-0066+L-0026+L-0027" + ("+L-0065" if bound.C_ub_sos_hi else ""),
        "domain": (
            "Pseudospectral T^3, N<=24, 2/3 dealias, Fourier support in shells "
            f"|k|^2 in {list(bound.support_cr0013 or [])}, single frozen polarization "
            f"(branch={bound.pol_branch} of _pol_basis) per wavevector. NOT continuum."
        ),
        "nu": 0.1,
        "energy": 0.5,
        "t_end": 0.02,
        "n_max": bound.n,
        "proposed_bound_M": bound.c0007_M,
        "shells": list(bound.support_cr0013 or []),
        "pol_branch": bound.pol_branch,
        "C_shor_sym": bound.C_shor_sym_cr0013,
        "C_ub_sos_hi": bound.C_ub_sos_hi,
        "C_dagger": bound.C_dagger,
        "D": bound.support_cr0013_D,
        "Omega_star": bound.Omega_star,
        "parent_open": "C-0007",
        "extends": "C-R-0008",
        "statement": (
            f"For divergence-free fields on N≤{bound.n} with Fourier support only on "
            f"shells |k|²∈{list(bound.support_cr0013 or [])} and a single frozen "
            f"polarization branch ({bound.pol_branch}) per k, dealiased Galerkin NS "
            f"(ν=0.1, E≤0.5) satisfies Ω(0.02)≤{bound.c0007_M}. Proof: one-pol Sym "
            f"Shor gives C≤{bound.C_shor_sym_cr0013}≤C_† ⇒ L-0027; high slab via L-0026."
            f"{sos_note} FINITE only."
        ),
        "clay_implication": (
            "None. One-pol structured subclass only; not all-IC C-0007; not Clay."
        ),
        "notes": (
            "Discovered by L-0064 one-pol scan; largest one-pol second-branch witness "
            f"under greedy growth (r≤20). Extends C-R-0008 by shells "
            f"{sorted(set(bound.support_cr0013 or []) - set(SUPPORT_CR0008))}."
        ),
    }


def save_lemma_l0066(
    bound: GalerkinBoundL0066,
    path: str | Path = "conjectures/proved_restricted/L-0066.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")


def save_cr0013(
    bound: GalerkinBoundL0066,
    path: str | Path = "conjectures/proved_restricted/C-R-0013.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cr0013_record(bound), indent=2), encoding="utf-8")
