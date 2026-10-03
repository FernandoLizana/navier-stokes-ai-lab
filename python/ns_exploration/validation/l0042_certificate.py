"""Certificate for L-0042 Sym-strengthening techo on {1..6}."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0040_sym_shor import sym_band_C_shor
from ns_exploration.conjectures.l0042_sym_strengthen_techo import (
    BAND_123456,
    ab_block_hybrid_C,
    column_triangle_C,
    lemma_l0042,
)


@dataclass
class CertificateL0042:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    C_sym_123456: float
    C_triangle_123456: float
    C_hybrid_AB: float
    C_rank1_lo_123456: float
    closes_123456: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0042_certificate(n: int = 24) -> CertificateL0042:
    b = lemma_l0042(n=n)
    return CertificateL0042(
        cert_id="CERT-L0042-sym-strengthen-techo-C0007-N24",
        lemma_id="L-0042",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        C_sym_123456=b.C_sym_123456,
        C_triangle_123456=b.C_triangle_123456,
        C_hybrid_AB=b.C_hybrid_AB,
        C_rank1_lo_123456=b.C_rank1_lo_123456,
        closes_123456=b.closes_123456,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0042_certificate(d: dict) -> tuple[bool, dict]:
    Cd = float(d["C_dagger"])
    c6 = sym_band_C_shor(int(d["n"]), BAND_123456)["C_shor_sym"]
    tri = column_triangle_C(int(d["n"]), BAND_123456)["C_triangle"]
    hy = ab_block_hybrid_C(int(d["n"]))["C_hybrid"]
    checks = {
        "sym_above_Cdagger": c6 > Cd,
        "sym_matches": abs(c6 - float(d["C_sym_123456"])) < 1e-6,
        "triangle_worse_than_sym": tri > c6,
        "hybrid_above_Cdagger": hy > Cd,
        "rank1_below_Cdagger": float(d["C_rank1_lo_123456"]) < Cd,
        "closes_false": d["closes_123456"] is False,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    return all(checks.values()), checks


def save_l0042_certificate(
    cert: CertificateL0042,
    path: str | Path = "certificates/CERT-L0042-sym-strengthen-techo-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
