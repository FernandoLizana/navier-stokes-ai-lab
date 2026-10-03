"""Certificate for L-0047 extended greedy physical fullsym + C-R-0012."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0047_greedy_fullsym import (
    SUPPORT_CR0012,
    lemma_l0047,
)
from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym


@dataclass
class CertificateL0047:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    support_cr0012: list[int]
    C_fullsym_cr0012: float
    closes_cr0012: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0047_certificate(
    n: int = 24, bound: object | None = None
) -> CertificateL0047:
    b = bound if bound is not None else lemma_l0047(n=n)
    return CertificateL0047(
        cert_id="CERT-L0047-greedy-fullsym-C0007-N24",
        lemma_id="L-0047",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        support_cr0012=list(b.support_cr0012 or []),
        C_fullsym_cr0012=b.C_fullsym_cr0012,
        closes_cr0012=b.closes_cr0012,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0047_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    n = int(d["n"])
    c = streaming_C_fullsym(n, SUPPORT_CR0012, dtype=np.float32)["C_fullsym"]
    checks = {
        "cr0012_le_Cdagger": c <= Cd + 1e-5,
        "cr0012_matches": abs(c - float(d["C_fullsym_cr0012"])) < 1e-4,
        "support_frozen": list(d["support_cr0012"]) == list(SUPPORT_CR0012),
        "includes_147": 147 in d["support_cr0012"],
        "n_shells_41": len(d["support_cr0012"]) == 41,
        "closes_cr": d["closes_cr0012"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0047_certificate(
    cert: CertificateL0047,
    path: str | Path = "certificates/CERT-L0047-greedy-fullsym-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
