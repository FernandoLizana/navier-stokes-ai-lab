"""Certificate for L-0041 shell-6 Sym + C-R-0007."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0040_sym_shor import sym_band_C_shor
from ns_exploration.conjectures.l0041_shell6_sym import (
    BAND_123456,
    SUPPORT_CR0007,
    lemma_l0041,
)


@dataclass
class CertificateL0041:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    band_123456_C: float
    seed_12346_C: float
    support_cr0007: list
    support_cr0007_C: float
    closes_cr0007: bool
    consecutive_6_techo: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0041_certificate(n: int = 24) -> CertificateL0041:
    b = lemma_l0041(n=n)
    return CertificateL0041(
        cert_id="CERT-L0041-shell6-sym-C0007-N24",
        lemma_id="L-0041",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        band_123456_C=b.band_123456_C,
        seed_12346_C=b.seed_12346_C,
        support_cr0007=list(b.support_cr0007 or SUPPORT_CR0007),
        support_cr0007_C=b.support_cr0007_C,
        closes_cr0007=b.closes_cr0007,
        consecutive_6_techo=b.consecutive_6_techo,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0041_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    c6 = sym_band_C_shor(int(d["n"]), BAND_123456)["C_shor_sym"]
    cr = sym_band_C_shor(int(d["n"]), d["support_cr0007"])["C_shor_sym"]
    checks = {
        "band6_above_Cdagger": c6 > Cd,
        "band6_matches": abs(c6 - float(d["band_123456_C"])) < 1e-6,
        "cr_le_Cdagger": cr <= Cd + 1e-6,
        "cr_matches": abs(cr - float(d["support_cr0007_C"])) < 1e-5,
        "seed_le_Cdagger": float(d["seed_12346_C"]) <= Cd,
        "cr0007_true": d["closes_cr0007"] is True,
        "techo_true": d["consecutive_6_techo"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0041_certificate(
    cert: CertificateL0041,
    path: str | Path = "certificates/CERT-L0041-shell6-sym-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
