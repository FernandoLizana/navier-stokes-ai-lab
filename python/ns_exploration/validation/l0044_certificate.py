"""Certificate for L-0044 physical Hermitian fullsym + C-R-0009."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0044_physical_fullsym import (
    SUPPORT_CR0009,
    lemma_l0044,
    physical_band_C_fullsym,
)


@dataclass
class CertificateL0044:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    support_cr0009: list[int]
    C_fullsym_123456: float
    C_sym_123456: float
    fft_tensor_err: float
    closes_cr0009: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0044_certificate(n: int = 24) -> CertificateL0044:
    b = lemma_l0044(n=n)
    return CertificateL0044(
        cert_id="CERT-L0044-physical-fullsym-C0007-N24",
        lemma_id="L-0044",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        support_cr0009=list(b.support_cr0009 or []),
        C_fullsym_123456=b.C_fullsym_123456,
        C_sym_123456=b.C_sym_123456,
        fft_tensor_err=b.fft_tensor_err,
        closes_cr0009=b.closes_cr0009,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0044_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    n = int(d["n"])
    band = physical_band_C_fullsym(n, SUPPORT_CR0009)
    checks = {
        "fullsym_le_Cdagger": band["C_fullsym"] <= Cd + 1e-9,
        "fullsym_matches": abs(band["C_fullsym"] - float(d["C_fullsym_123456"])) < 1e-6,
        "sym_above_Cdagger": band["C_sym"] > Cd,
        "support_frozen": list(d["support_cr0009"]) == list(SUPPORT_CR0009),
        "fft_err_small": float(d["fft_tensor_err"]) < 1e-11,
        "closes_cr": d["closes_cr0009"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0044_certificate(
    cert: CertificateL0044,
    path: str | Path = "certificates/CERT-L0044-physical-fullsym-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
