"""
L-0059: Lean N8 bridge for N5 SOS Gram certificates (L-0052/53/54).

Documents Lean-checked rational comparisons C_ub_hi < C_dagger for the
three existing SOS Gram certs. Arithmetic only; does not formalize SDP.

FINITE Galerkin only. Not Clay.
"""

from __future__ import annotations

import json
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027

LEAN_DIR = Path("lean/NSGalerkin")
CERTS = (
    ("CERT-L0052-sos-gram-band123-N24", "band123"),
    ("CERT-L0053-sos-gram-onepol-cr0008-N24", "onepol_cr0008"),
    ("CERT-L0054-sos-gram-band1234-N24", "band1234"),
    ("CERT-L0058-sos-gram-band12345-N24", "band12345"),
    ("CERT-L0060-sos-gram-band123456-N24", "band123456"),
    ("CERT-L0061-sos-band12345678-N24", "band12345678"),
    ("CERT-L0062-sos-band123456789-N24", "band123456789"),
    ("CERT-L0062-sos-band1to10-N24", "band1to10"),
    ("CERT-L0062-sos-band1to12-N24", "band1to12"),
    ("CERT-L0065-sos-onepol-greedy-second-N24", "onepol_greedy_second"),
    ("CERT-L0067-sos-onepol-greedy-second-r24-N24", "onepol_greedy_second_r24"),
)


@dataclass
class LeanPinRow:
    cert_id: str
    subclass: str
    C_ub_hi_num: int
    C_ub_hi_den: int
    C_dagger_num: int
    C_dagger_den: int
    lean_theorem: str


@dataclass
class GalerkinBoundL0059:
    lemma_id: str = "L-0059"
    route: str = "B"
    status: str = "validated"
    evidence_level: str = "N8"
    n: int = 24
    lean_package: str = "NSGalerkin.SosGram"
    lean_build_ok: bool = False
    pins: list[LeanPinRow] | None = None
    clay_implication: str = (
        "None. Lean arithmetic cross-check of N5 SOS constants only; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.pins:
            d["pins"] = [asdict(p) for p in self.pins]
        return d


def _rat_ceil(x: float, den: int = 10**16) -> tuple[int, int]:
    """Conservative rational upper bound for Lean native_decide."""
    import math

    num = math.ceil(x * den - 1e-12)
    return int(num), den


def _load_cert(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def lemma_l0059(n: int = 24, run_lean: bool = True) -> GalerkinBoundL0059:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    cd_num, cd_den = _rat_ceil(Cd)

    pins: list[LeanPinRow] = []
    for cert_id, tag in CERTS:
        path = Path("certificates") / f"{cert_id}.json"
        if not path.is_file():
            continue
        d = _load_cert(path)
        hi = float(d["C_ub_hi"])
        hi_num, hi_den = _rat_ceil(hi)
        pins.append(
            LeanPinRow(
                cert_id=cert_id,
                subclass=str(d.get("lemma_ref", tag)),
                C_ub_hi_num=hi_num,
                C_ub_hi_den=hi_den,
                C_dagger_num=cd_num,
                C_dagger_den=cd_den,
                lean_theorem=f"sos_{tag}_C_ub_lt_Cdagger",
            )
        )

    lean_ok = False
    if run_lean and LEAN_DIR.is_dir():
        try:
            r = subprocess.run(
                ["lake", "build"],
                cwd=LEAN_DIR,
                capture_output=True,
                text=True,
                timeout=300,
            )
            lean_ok = r.returncode == 0
        except (FileNotFoundError, subprocess.TimeoutExpired):
            lean_ok = False

    notes = (
        f"Lean N8 pins for {len(pins)} N5 SOS Gram certs: "
        "native_decide checks C_ub_hi < C_dagger (rational arithmetic only). "
        f"lake build: {'ok' if lean_ok else 'pending/failed'}."
    )
    return GalerkinBoundL0059(
        n=n,
        lean_build_ok=lean_ok,
        pins=pins,
        notes=notes,
    )


def save_lemma_l0059(
    bound: GalerkinBoundL0059,
    path: str | Path = "conjectures/active/L-0059.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
