"""Certificate for L-0039 sparse support Shor + C-R-0005."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0038_band_svd_shor import exact_band_C_shor
from ns_exploration.conjectures.l0039_sparse_support_shor import (
    QUAD_1346,
    SUPPORT_A,
    lemma_l0039,
    sparse_band_C_shor,
)


@dataclass
class CertificateL0039:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    quad_1346_C: float
    support_A: list
    support_A_C: float
    support_B_C: float
    closes_cr0005: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0039_certificate(n: int = 24) -> CertificateL0039:
    b = lemma_l0039(n=n)
    return CertificateL0039(
        cert_id="CERT-L0039-sparse-support-C0007-N24",
        lemma_id="L-0039",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        quad_1346_C=b.quad_1346_C,
        support_A=list(b.support_A or SUPPORT_A),
        support_A_C=b.support_A_C,
        support_B_C=b.support_B_C,
        closes_cr0005=b.closes_cr0005,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0039_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    q = exact_band_C_shor(int(d["n"]), QUAD_1346)["C_shor"]
    a = sparse_band_C_shor(int(d["n"]), d["support_A"])["C_shor"]
    checks = {
        "quad_le_Cdagger": q <= Cd + 1e-8,
        "quad_matches": abs(q - float(d["quad_1346_C"])) < 1e-6,
        "supportA_le_Cdagger": a <= Cd + 1e-5,
        "supportA_matches": abs(a - float(d["support_A_C"])) < 1e-4,
        "supportB_le_Cdagger": float(d["support_B_C"]) <= Cd + 1e-5,
        "cr0005_true": d["closes_cr0005"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0039_certificate(
    cert: CertificateL0039,
    path: str | Path = "certificates/CERT-L0039-sparse-support-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
