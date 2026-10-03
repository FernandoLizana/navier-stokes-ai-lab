"""Certificate for L-0036 signed Shor stretch techo."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0036_signed_shor_stretch import (
    lemma_l0036,
    signed_flattening_C_lower,
)


@dataclass
class CertificateL0036:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    D: int
    n_pairs: int
    L_op_lower: float
    C_shor_lower: float
    closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0036_certificate(n: int = 12, iters: int = 10) -> CertificateL0036:
    b = lemma_l0036(n=n, iters=iters, include_n24_probe=False)
    return CertificateL0036(
        cert_id="CERT-L0036-signed-shor-C0007-N12",
        lemma_id="L-0036",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        D=b.D,
        n_pairs=b.n_pairs,
        L_op_lower=b.L_op_lower,
        C_shor_lower=b.C_shor_lower,
        closes_c0007=b.closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0036_certificate(d: dict) -> tuple[bool, dict]:
    probe = signed_flattening_C_lower(n=int(d["n"]), iters=8)
    checks = {
        "D_match": probe["D"] == int(d["D"]),
        "pairs_match": probe["n_pairs"] == int(d["n_pairs"]),
        "C_shor_above_Cdagger": float(d["C_shor_lower"]) > float(d["C_dagger"]),
        "recomputed_also_above": float(probe["C_shor_lower"]) > float(d["C_dagger"]),
        "closes_false": d["closes_c0007"] is False,
    }
    return all(checks.values()), checks


def save_l0036_certificate(
    cert: CertificateL0036,
    path: str | Path = "certificates/CERT-L0036-signed-shor-C0007-N12.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
