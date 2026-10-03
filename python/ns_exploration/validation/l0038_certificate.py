"""Certificate for L-0038 band SVD Shor + C-R-0004."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0038_band_svd_shor import (
    exact_band_C_shor,
    lemma_l0038,
)


@dataclass
class CertificateL0038:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    band_123_C: float
    band_1234_C: float
    band_12_C: float
    band_125_C: float
    closes_cr0004: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0038_certificate(n: int = 24) -> CertificateL0038:
    b = lemma_l0038(n=n)
    return CertificateL0038(
        cert_id="CERT-L0038-band-svd-C0007-N24",
        lemma_id="L-0038",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        band_123_C=b.band_123_C,
        band_1234_C=b.band_1234_C,
        band_12_C=b.band_12_C,
        band_125_C=b.band_125_C,
        closes_cr0004=b.closes_cr0004,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0038_certificate(d: dict) -> tuple[bool, dict]:
    c123 = exact_band_C_shor(int(d["n"]), (1, 2, 3))["C_shor"]
    c1234 = exact_band_C_shor(int(d["n"]), (1, 2, 3, 4))["C_shor"]
    Cd = float(d["C_dagger"])
    checks = {
        "band123_matches": abs(c123 - float(d["band_123_C"])) < 1e-6,
        "band123_le_Cdagger": c123 <= Cd + 1e-9,
        "band1234_above_Cdagger": c1234 > Cd,
        "band125_le_Cdagger": float(d["band_125_C"]) <= Cd + 1e-9,
        "cr0004_true": d["closes_cr0004"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0038_certificate(
    cert: CertificateL0038,
    path: str | Path = "certificates/CERT-L0038-band-svd-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
