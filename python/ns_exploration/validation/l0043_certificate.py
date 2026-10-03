"""Certificate for L-0043 one-pol Sym + shell-block/Ky-Fan techos."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0043_onepol_shellblock import (
    BAND_ONEPOL_TECHO,
    SUPPORT_CR0008,
    kyfan_partial_C,
    lemma_l0043,
    shellblock_C_lower,
    sym_onepol_C_shor,
)


@dataclass
class CertificateL0043:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    support_cr0008: list[int]
    support_cr0008_C: float
    onepol_techo_12345689_C: float
    C_shellblock_lo: float
    C_kyfan_partial: float
    closes_cr0008: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0043_certificate(n: int = 24) -> CertificateL0043:
    b = lemma_l0043(n=n)
    return CertificateL0043(
        cert_id="CERT-L0043-onepol-shellblock-C0007-N24",
        lemma_id="L-0043",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        support_cr0008=list(b.support_cr0008 or []),
        support_cr0008_C=b.support_cr0008_C,
        onepol_techo_12345689_C=b.onepol_techo_12345689_C,
        C_shellblock_lo=b.C_shellblock_lo,
        C_kyfan_partial=b.C_kyfan_partial,
        closes_cr0008=b.closes_cr0008,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0043_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    n = int(d["n"])
    cr = sym_onepol_C_shor(n, SUPPORT_CR0008)["C_shor_sym"]
    techo = sym_onepol_C_shor(n, BAND_ONEPOL_TECHO)["C_shor_sym"]
    sb = shellblock_C_lower(n=n)["C_shellblock_lo"]
    kf = kyfan_partial_C(n=n, n_terms=4)["C_kyfan_partial"]
    checks = {
        "cr0008_le_Cdagger": cr <= Cd + 1e-9,
        "cr0008_matches": abs(cr - float(d["support_cr0008_C"])) < 1e-6,
        "support_frozen": list(d["support_cr0008"]) == list(SUPPORT_CR0008),
        "onepol_techo_above": techo > Cd,
        "techo_matches": abs(techo - float(d["onepol_techo_12345689_C"])) < 1e-5,
        "shellblock_above": sb > Cd,
        "kyfan_above": kf > Cd,
        "closes_cr": d["closes_cr0008"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0043_certificate(
    cert: CertificateL0043,
    path: str | Path = "certificates/CERT-L0043-onepol-shellblock-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
