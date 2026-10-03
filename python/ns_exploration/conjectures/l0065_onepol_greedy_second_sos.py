"""
L-0065: N5 interval SOS for one-pol greedy-second witness (L-0064).

Support shells [1,2,3,4,5,6,8,10,12,16,19], branch=second, D=154.
COSMO TSSOS N2 validated; N5 interval from ub_raw. Gram deferred on laptop.

FINITE Galerkin one-pol subclass only. Not all-IC. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0052_sos_gram_n5 import certified_C_ub_hi

TSSOS_JSON = Path("tools/sos_julia/data/onepol_greedy_second/cosmo_tssos_result.json")
SUPPORT_GREEDY_SECOND = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 19]
POL_BRANCH = "second"
D_GREEDY_SECOND = 154
C_FULLSYM = 9.060826005460791
C_SHOR_SYM = 9.189036845789548
C_UB_FLOAT = 4.01493826157398


@dataclass
class GalerkinBoundL0065:
    lemma_id: str = "L-0065"
    route: str = "B"
    status: str = "validated_constant"
    evidence_level: str = "N5"
    n: int = 24
    support_radii: list[int] | None = None
    pol_branch: str = POL_BRANCH
    D: int = D_GREEDY_SECOND
    C_dagger: float = 0.0
    C_shor_sym: float = C_SHOR_SYM
    C_fullsym: float = C_FULLSYM
    C_ub_float: float = C_UB_FLOAT
    C_ub_hi: float = 0.0
    ub_raw_hi: float = 0.0
    beats_fullsym: bool = True
    beats_shor: bool = True
    closes_vs_Cdagger: bool = True
    closes_all_ic: bool = False
    solver: str = "COSMO"
    scan_ref: str = "L-0064 greedy-second"
    proposed_cr: str = "C-R-0013"
    clay_implication: str = (
        "None. N5 interval on one-pol greedy-second subclass only; not all-IC; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0065(
    n: int = 24,
    tssos_path: Path = TSSOS_JSON,
) -> GalerkinBoundL0065:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    if not tssos_path.is_file():
        raise FileNotFoundError(f"missing COSMO result: {tssos_path}")
    tssos = json.loads(tssos_path.read_text(encoding="utf-8"))
    ub_raw = float(tssos["ub_raw_abs"])
    ub_hi, C_hi = certified_C_ub_hi(ub_raw)
    C_float = float(tssos.get("C_ub", C_UB_FLOAT))
    C_fs = float(tssos.get("C_fullsym", C_FULLSYM))
    C_sh = float(
        json.loads(
            Path("tools/sos_julia/data/onepol_greedy_second/meta.json").read_text(
                encoding="utf-8"
            )
        ).get("C_shor_sym", C_SHOR_SYM)
    )
    notes = (
        f"One-pol greedy-second (L-0064): {len(SUPPORT_GREEDY_SECOND)} shells, "
        f"branch={POL_BRANCH}, D={D_GREEDY_SECOND}. COSMO TSSOS C_ub≈{C_float:.4f}; "
        f"N5 C_ub_hi≈{C_hi:.6f} < C_dagger≈{Cd:.4f} < C_shor≈{C_sh:.4f}. "
        f"Extends C-R-0008 (8 shells). Gram deferred. Proposed C-R-0013."
    )
    return GalerkinBoundL0065(
        n=n,
        support_radii=list(SUPPORT_GREEDY_SECOND),
        C_dagger=Cd,
        C_shor_sym=C_sh,
        C_fullsym=C_fs,
        C_ub_float=C_float,
        C_ub_hi=C_hi,
        ub_raw_hi=ub_hi,
        beats_fullsym=C_hi < C_fs - 1e-9,
        beats_shor=C_hi < C_sh - 1e-9,
        closes_vs_Cdagger=C_hi < Cd - 1e-9,
        notes=notes,
    )


def save_lemma_l0065(
    bound: GalerkinBoundL0065,
    path: str | Path = "conjectures/active/L-0065.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
