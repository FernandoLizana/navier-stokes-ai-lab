"""Certificate for L-0066 / C-R-0013 one-pol greedy-second."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0043_onepol_shellblock import SUPPORT_CR0008, sym_onepol_C_shor
from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import SUPPORT_GREEDY_SECOND, POL_BRANCH
from ns_exploration.conjectures.l0066_cr0013_onepol_greedy import lemma_l0066


@dataclass
class CertificateL0066:
    cert_id: str
    lemma_id: str
    restricted_id: str
    n: int
    c0007_M: float
    C_dagger: float
    support_cr0013: list[int]
    pol_branch: str
    C_shor_sym: float
    C_ub_sos_hi: float | None
    closes_cr0013: bool
    extends_cr0008: bool
    closes_c0007_all_ic: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0066_certificate(n: int = 24) -> CertificateL0066:
    b = lemma_l0066(n=n)
    return CertificateL0066(
        cert_id="CERT-L0066-onepol-greedy-second-C0007-N24",
        lemma_id="L-0066",
        restricted_id="C-R-0013",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        support_cr0013=list(b.support_cr0013 or []),
        pol_branch=b.pol_branch,
        C_shor_sym=b.C_shor_sym_cr0013,
        C_ub_sos_hi=b.C_ub_sos_hi,
        closes_cr0013=b.closes_cr0013,
        extends_cr0008=b.extends_cr0008,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        clay_implication=b.clay_implication,
    )


def verify_l0066_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    n = int(d["n"])
    Cd = float(d["C_dagger"])
    ref = sym_onepol_C_shor(n, SUPPORT_GREEDY_SECOND, branch=POL_BRANCH)
    C = float(ref["C_shor_sym"])
    checks = {
        "support_frozen": list(d["support_cr0013"]) == list(SUPPORT_GREEDY_SECOND),
        "C_le_Cdagger": C <= Cd + 1e-9,
        "C_matches": abs(C - float(d["C_shor_sym"])) < 1e-6,
        "extends_cr0008": set(SUPPORT_CR0008).issubset(set(SUPPORT_GREEDY_SECOND)),
        "closes_cr": d["closes_cr0013"] is True,
        "all_ic_false": d["closes_c0007_all_ic"] is False,
    }
    if d.get("C_ub_sos_hi") is not None:
        checks["sos_below_Cdagger"] = float(d["C_ub_sos_hi"]) < Cd - 1e-9
    return all(checks.values()), checks


def save_l0066_certificate(
    cert: CertificateL0066,
    path: str | Path = "certificates/CERT-L0066-onepol-greedy-second-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
