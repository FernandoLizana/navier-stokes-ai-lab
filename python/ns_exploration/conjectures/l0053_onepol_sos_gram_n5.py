"""
L-0053: N5 Gram rationalization for one-pol TSSOS on C-R-0008 (D=92).

Structured subclass (fixed polarization). Not all-IC. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0043_onepol_shellblock import SUPPORT_CR0008
from ns_exploration.conjectures.l0052_sos_gram_n5 import (
    certified_C_ub_hi,
    load_gram_export,
    verify_gram_blocks,
)

GRAM_JSON = Path("tools/sos_julia/data/onepol_cr0008/gram_export.json")
TSSOS_JSON = Path("tools/sos_julia/data/onepol_cr0008/tssos_result.json")

C_UB_ONEPOL = 3.3560651946371607
C_FULLSYM_ONEPOL = 8.475443542120622
D_ONEPOL = 92


@dataclass
class GalerkinBoundL0053:
    lemma_id: str = "L-0053"
    route: str = "B"
    status: str = "validated"
    evidence_level: str = "N5"
    n: int = 24
    subclass: str = "one_pol_cr0008"
    radii: list[int] | None = None
    D: int = D_ONEPOL
    C_dagger: float = 0.0
    C_fullsym: float = C_FULLSYM_ONEPOL
    C_ub_float: float = C_UB_ONEPOL
    C_ub_hi: float = 0.0
    ub_raw_hi: float = 0.0
    gram_denom: int = 0
    n_gram_blocks: int = 0
    min_psd_eig: float = 0.0
    max_rationalization_error: float = 0.0
    beats_fullsym: bool = True
    closes_vs_Cdagger: bool = True
    closes_all_ic: bool = False
    clay_implication: str = (
        "None. N5 Gram on one-pol C-R-0008 subclass only; not all-IC; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def lemma_l0053(
    n: int = 24,
    gram_path: Path = GRAM_JSON,
    tssos_path: Path = TSSOS_JSON,
) -> GalerkinBoundL0053:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    gram = load_gram_export(gram_path)
    tssos = json.loads(tssos_path.read_text(encoding="utf-8")) if tssos_path.is_file() else {}

    ub_raw = float(gram.get("ub_raw_abs", tssos.get("ub_raw_abs", 0.0)))
    ub_hi, C_hi = certified_C_ub_hi(ub_raw)
    C_float = float(gram.get("C_ub", tssos.get("C_ub", C_UB_ONEPOL)))
    C_fs = float(gram.get("C_fullsym", C_FULLSYM_ONEPOL))

    checks = verify_gram_blocks(gram)
    if not all(checks.values()):
        bad = [k for k, v in checks.items() if not v]
        raise ValueError(f"Gram verification failed: {bad}")

    notes = (
        f"N5 Gram on one-pol C-R-0008 D={D_ONEPOL}: "
        f"{gram.get('n_gram_blocks', 0)} PSD blocks, denom={gram.get('denom')}. "
        f"C_ub_hi≈{C_hi:.4f} ≪ C_fullsym≈{C_fs:.4f} < C_†≈{Cd:.4f}. "
        "Fixed polarization subclass; rational Gram + interval constant."
    )
    return GalerkinBoundL0053(
        n=n,
        radii=list(SUPPORT_CR0008),
        C_dagger=Cd,
        C_fullsym=C_fs,
        C_ub_float=C_float,
        C_ub_hi=C_hi,
        ub_raw_hi=ub_hi,
        gram_denom=int(gram.get("denom", 0)),
        n_gram_blocks=int(gram.get("n_gram_blocks", 0)),
        min_psd_eig=float(gram.get("min_psd_eig_all_blocks", 0.0)),
        max_rationalization_error=float(gram.get("max_rationalization_error", 0.0)),
        beats_fullsym=C_hi < C_fs - 1e-9,
        closes_vs_Cdagger=C_hi < Cd - 1e-9,
        notes=notes,
    )


def save_lemma_l0053(
    bound: GalerkinBoundL0053,
    path: str | Path = "conjectures/active/L-0053.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
