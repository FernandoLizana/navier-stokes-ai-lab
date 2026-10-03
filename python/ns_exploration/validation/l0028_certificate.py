"""Certificate for L-0028 embedding techo."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0002_viscous_galerkin import stretch_constant
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0028_embedding_techo import (
    M_star_for_Cdagger,
    lemma_l0028,
    modes_upto_radius,
)


@dataclass
class CertificateL0028:
    cert_id: str
    lemma_id: str
    n: int
    C_dagger: float
    M_star: float
    M_full: int
    C_full_L0002: float
    M_shell_r1: int
    C_shell_r1: float
    M_hard_R2: int
    C_hard_R2: float
    embedding_can_meet_Cdagger: bool
    proved_closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0028_certificate(n: int = 24) -> CertificateL0028:
    b = lemma_l0028(n=n, empirical=False)
    return CertificateL0028(
        cert_id="CERT-L0028-embedding-techo-C0007-N24",
        lemma_id="L-0028",
        n=b.n,
        C_dagger=b.C_dagger,
        M_star=b.M_star,
        M_full=b.M_full,
        C_full_L0002=b.C_full_L0002,
        M_shell_r1=b.M_shell_r1,
        C_shell_r1=b.C_shell_r1,
        M_hard_R2=b.M_hard_R2,
        C_hard_R2=b.C_hard_R2,
        embedding_can_meet_Cdagger=b.embedding_can_meet_Cdagger,
        proved_closes_c0007=b.proved_closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0028_certificate(d: dict) -> tuple[bool, dict]:
    n = int(d["n"])
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    Mstar = M_star_for_Cdagger(Cd)
    M_full = modes_upto_radius(n, 10**9)
    M_r1 = modes_upto_radius(n, 1)
    M_R2 = modes_upto_radius(n, 2)
    checks = {
        "C_dagger_matches": abs(Cd - float(d["C_dagger"])) < 1e-6,
        "M_star_matches": abs(Mstar - float(d["M_star"])) < 1e-6,
        "M_full_matches": M_full == int(d["M_full"]),
        "C_full_matches": abs(stretch_constant(M_full) - float(d["C_full_L0002"]))
        < 1e-6,
        "M_r1_is_6": M_r1 == 6 and int(d["M_shell_r1"]) == 6,
        "C_r1_above_Cdagger": stretch_constant(M_r1) > Cd,
        "C_R2_above_Cdagger": stretch_constant(M_R2) > Cd,
        "M_star_below_6": Mstar < 6,
        "embedding_cannot_meet": d["embedding_can_meet_Cdagger"] is False,
        "proved_closes_false": d["proved_closes_c0007"] is False,
    }
    return all(checks.values()), checks


def save_l0028_certificate(
    cert: CertificateL0028,
    path: str | Path = "certificates/CERT-L0028-embedding-techo-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
