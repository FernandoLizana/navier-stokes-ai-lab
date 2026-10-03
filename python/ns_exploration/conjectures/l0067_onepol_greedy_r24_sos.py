"""
L-0067: N5 interval SOS for one-pol greedy-second + shell 24 (C-R-0014).

12 shells, branch=second, D=178. Extends C-R-0013 by shell 24 (L-0064 scan).

FINITE Galerkin one-pol subclass only. Not all-IC. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0052_sos_gram_n5 import certified_C_ub_hi
from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import SUPPORT_GREEDY_SECOND

SUPPORT_CR0014 = SUPPORT_GREEDY_SECOND + [24]
POL_BRANCH = "second"
D_CR0014 = 178
TSSOS_JSON = Path("tools/sos_julia/data/onepol_greedy_second_r24/cosmo_tssos_result.json")
META_JSON = Path("tools/sos_julia/data/onepol_greedy_second_r24/meta.json")


@dataclass
class GalerkinBoundL0067:
    lemma_id: str = "L-0067"
    route: str = "B"
    status: str = "validated_constant"
    evidence_level: str = "N5"
    n: int = 24
    support_radii: list[int] | None = None
    pol_branch: str = POL_BRANCH
    D: int = D_CR0014
    C_dagger: float = 0.0
    C_shor_sym: float = 0.0
    C_fullsym: float = 0.0
    C_ub_float: float = 0.0
    C_ub_hi: float = 0.0
    ub_raw_hi: float = 0.0
    beats_shor: bool = True
    closes_vs_Cdagger: bool = True
    closes_all_ic: bool = False
    extends_cr: str = "C-R-0013"
    proposed_cr: str = "C-R-0014"
    solver: str = "COSMO"
    clay_implication: str = (
        "None. N5 interval on one-pol 12-shell witness; not all-IC; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0067(n: int = 24, tssos_path: Path = TSSOS_JSON) -> GalerkinBoundL0067:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    if not tssos_path.is_file():
        raise FileNotFoundError(f"missing COSMO result: {tssos_path}")
    tssos = json.loads(tssos_path.read_text(encoding="utf-8"))
    meta = json.loads(META_JSON.read_text(encoding="utf-8")) if META_JSON.is_file() else {}
    ub_raw = float(tssos["ub_raw_abs"])
    ub_hi, C_hi = certified_C_ub_hi(ub_raw)
    C_float = float(tssos.get("C_ub", 0.0))
    C_fs = float(tssos.get("C_fullsym", meta.get("C_fullsym", 0.0)))
    C_sh = float(meta.get("C_shor_sym", 0.0))
    notes = (
        f"One-pol greedy-second+24: {len(SUPPORT_CR0014)} shells D={D_CR0014}. "
        f"COSMO C_ub≈{C_float:.4f}; N5 C_ub_hi≈{C_hi:.6f} < C_dagger≈{Cd:.4f} "
        f"< C_shor≈{C_sh:.4f}. Extends C-R-0013. Proposed C-R-0014."
    )
    return GalerkinBoundL0067(
        n=n,
        support_radii=list(SUPPORT_CR0014),
        C_dagger=Cd,
        C_shor_sym=C_sh,
        C_fullsym=C_fs,
        C_ub_float=C_float,
        C_ub_hi=C_hi,
        ub_raw_hi=ub_hi,
        beats_shor=C_hi < C_sh - 1e-9,
        closes_vs_Cdagger=C_hi < Cd - 1e-9,
        notes=notes,
    )


def save_lemma_l0067(
    bound: GalerkinBoundL0067,
    path: str | Path = "conjectures/active/L-0067.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
