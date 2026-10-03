"""Certificate for L-0040 Sym Shor + C-R-0006."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0040_sym_shor import (
    SUPPORT_CR0006,
    lemma_l0040,
    sym_band_C_shor,
)


@dataclass
class CertificateL0040:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    band_12345_C: float
    band_123456_C: float
    band_1234_C_sym: float
    greedy9_C: float
    closes_cr0006: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0040_certificate(n: int = 24) -> CertificateL0040:
    b = lemma_l0040(n=n)
    return CertificateL0040(
        cert_id="CERT-L0040-sym-shor-C0007-N24",
        lemma_id="L-0040",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        band_12345_C=b.band_12345_C,
        band_123456_C=b.band_123456_C,
        band_1234_C_sym=b.band_1234_C_sym,
        greedy9_C=b.greedy9_C,
        closes_cr0006=b.closes_cr0006,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0040_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    c5 = sym_band_C_shor(int(d["n"]), SUPPORT_CR0006)["C_shor_sym"]
    c6 = sym_band_C_shor(int(d["n"]), (1, 2, 3, 4, 5, 6))["C_shor_sym"]
    checks = {
        "band5_matches": abs(c5 - float(d["band_12345_C"])) < 1e-6,
        "band5_le_Cdagger": c5 <= Cd + 1e-9,
        "band6_above_Cdagger": c6 > Cd,
        "band4_le_Cdagger": float(d["band_1234_C_sym"]) <= Cd,
        "greedy9_le_Cdagger": float(d["greedy9_C"]) <= Cd + 1e-6,
        "cr0006_true": d["closes_cr0006"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0040_certificate(
    cert: CertificateL0040,
    path: str | Path = "certificates/CERT-L0040-sym-shor-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
