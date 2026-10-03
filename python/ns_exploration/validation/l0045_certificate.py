"""Certificate for L-0045 streaming physical fullsym + C-R-0010."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0045_physical_band import (
    SUPPORT_CR0010,
    lemma_l0045,
    radii_upto,
)
from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym


@dataclass
class CertificateL0045:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    support_cr0010: list[int]
    C_fullsym_rmax25: float
    C_fullsym_rmax26: float
    closes_cr0010: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0045_certificate(
    n: int = 24, bound: object | None = None
) -> CertificateL0045:
    b = bound if bound is not None else lemma_l0045(n=n)
    return CertificateL0045(
        cert_id="CERT-L0045-streaming-fullsym-C0007-N24",
        lemma_id="L-0045",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        support_cr0010=list(b.support_cr0010 or []),
        C_fullsym_rmax25=b.C_fullsym_rmax25,
        C_fullsym_rmax26=b.C_fullsym_rmax26,
        closes_cr0010=b.closes_cr0010,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0045_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    n = int(d["n"])
    c25 = streaming_C_fullsym(n, SUPPORT_CR0010)["C_fullsym"]
    c26 = streaming_C_fullsym(n, radii_upto(n, 26))["C_fullsym"]
    checks = {
        "r25_le_Cdagger": c25 <= Cd + 1e-9,
        "r25_matches": abs(c25 - float(d["C_fullsym_rmax25"])) < 1e-5,
        "r26_above": c26 > Cd,
        "r26_matches": abs(c26 - float(d["C_fullsym_rmax26"])) < 1e-5,
        "support_frozen": list(d["support_cr0010"]) == list(SUPPORT_CR0010),
        "closes_cr": d["closes_cr0010"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0045_certificate(
    cert: CertificateL0045,
    path: str | Path = "certificates/CERT-L0045-streaming-fullsym-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
