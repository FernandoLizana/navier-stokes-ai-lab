"""Certificate for L-0029 Frobenius stretch embedding."""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0002_viscous_galerkin import stretch_constant
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0028_embedding_techo import modes_upto_radius
from ns_exploration.conjectures.l0029_frobenius_stretch import (
    M_star_for_Cdagger_frobenius,
    lemma_l0029,
    stretch_constant_frobenius,
)


@dataclass
class CertificateL0029:
    cert_id: str
    lemma_id: str
    n: int
    C_dagger: float
    M_full: int
    C_F_full: float
    C_L0002_full: float
    M_star_F: float
    M_shell_r1: int
    C_F_shell_r1: float
    shell_r1_meets_Cdagger: bool
    shell_r1_reaches_low_slab_edge: bool
    frobenius_closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0029_certificate(n: int = 24) -> CertificateL0029:
    b = lemma_l0029(n=n)
    return CertificateL0029(
        cert_id="CERT-L0029-frobenius-stretch-C0007-N24",
        lemma_id="L-0029",
        n=b.n,
        C_dagger=b.C_dagger,
        M_full=b.M_full,
        C_F_full=b.C_F_full,
        C_L0002_full=b.C_L0002_full,
        M_star_F=b.M_star_F,
        M_shell_r1=b.M_shell_r1,
        C_F_shell_r1=b.C_F_shell_r1,
        shell_r1_meets_Cdagger=b.shell_r1_meets_Cdagger,
        shell_r1_reaches_low_slab_edge=b.shell_r1_reaches_low_slab_edge,
        frobenius_closes_c0007=b.frobenius_closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0029_certificate(d: dict) -> tuple[bool, dict]:
    n = int(d["n"])
    b7 = lemma_l0027(n=n, empirical=False)
    Cd = b7.C_dagger
    M_full = modes_upto_radius(n, 10**9)
    M_r1 = modes_upto_radius(n, 1)
    C_F = stretch_constant_frobenius(M_full)
    C_old = stretch_constant(M_full)
    Mstar = M_star_for_Cdagger_frobenius(Cd)
    C_r1 = stretch_constant_frobenius(M_r1)
    checks = {
        "C_dagger_matches": abs(Cd - float(d["C_dagger"])) < 1e-6,
        "M_full_matches": M_full == int(d["M_full"]),
        "C_F_matches": abs(C_F - float(d["C_F_full"])) < 1e-6,
        "C_L0002_matches": abs(C_old - float(d["C_L0002_full"])) < 1e-6,
        "factor_sqrt3": abs(C_old / C_F - math.sqrt(3.0)) < 1e-9,
        "M_star_F_matches": abs(Mstar - float(d["M_star_F"])) < 1e-6,
        "C_r1_le_Cdagger": C_r1 <= Cd + 1e-12,
        "M_star_F_lt_12": Mstar < 12,
        "C_F_full_gt_Cdagger": C_F > Cd,
        "r1_no_edge": d["shell_r1_reaches_low_slab_edge"] is False,
        "closes_false": d["frobenius_closes_c0007"] is False,
    }
    return all(checks.values()), checks


def save_l0029_certificate(
    cert: CertificateL0029,
    path: str | Path = "certificates/CERT-L0029-frobenius-stretch-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
