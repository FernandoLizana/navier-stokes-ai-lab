"""Certificate for L-0046 greedy physical fullsym + C-R-0011."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0046_greedy_fullsym import (
    SUPPORT_CR0011,
    lemma_l0046,
)
from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym


@dataclass
class CertificateL0046:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    support_cr0011: list[int]
    C_fullsym_cr0011: float
    closes_cr0011: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0046_certificate(
    n: int = 24, bound: object | None = None
) -> CertificateL0046:
    b = bound if bound is not None else lemma_l0046(n=n)
    return CertificateL0046(
        cert_id="CERT-L0046-greedy-fullsym-C0007-N24",
        lemma_id="L-0046",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        support_cr0011=list(b.support_cr0011 or []),
        C_fullsym_cr0011=b.C_fullsym_cr0011,
        closes_cr0011=b.closes_cr0011,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0046_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    n = int(d["n"])
    c = streaming_C_fullsym(n, SUPPORT_CR0011)["C_fullsym"]
    checks = {
        "cr0011_le_Cdagger": c <= Cd + 1e-9,
        "cr0011_matches": abs(c - float(d["C_fullsym_cr0011"])) < 1e-5,
        "support_frozen": list(d["support_cr0011"]) == list(SUPPORT_CR0011),
        "includes_147": 147 in d["support_cr0011"],
        "closes_cr": d["closes_cr0011"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0046_certificate(
    cert: CertificateL0046,
    path: str | Path = "certificates/CERT-L0046-greedy-fullsym-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
