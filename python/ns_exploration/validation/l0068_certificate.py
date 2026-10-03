"""Certificate for L-0068 / C-R-0014."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0067_onepol_greedy_r24_sos import SUPPORT_CR0014, POL_BRANCH
from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import SUPPORT_GREEDY_SECOND
from ns_exploration.conjectures.l0043_onepol_shellblock import sym_onepol_C_shor
from ns_exploration.conjectures.l0068_cr0014_onepol_greedy import lemma_l0068


@dataclass
class CertificateL0068:
    cert_id: str
    lemma_id: str
    restricted_id: str
    n: int
    support_cr0014: list[int]
    C_shor_sym: float
    C_ub_sos_hi: float | None
    closes_cr0014: bool
    extends_cr0013: bool
    closes_c0007_all_ic: bool

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0068_certificate(n: int = 24) -> CertificateL0068:
    b = lemma_l0068(n=n)
    return CertificateL0068(
        cert_id="CERT-L0068-onepol-greedy-second-r24-C0007-N24",
        lemma_id="L-0068",
        restricted_id="C-R-0014",
        n=n,
        support_cr0014=list(b.support_cr0014 or []),
        C_shor_sym=b.C_shor_sym_cr0014,
        C_ub_sos_hi=b.C_ub_sos_hi,
        closes_cr0014=b.closes_cr0014,
        extends_cr0013=b.extends_cr0013,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
    )


def verify_l0068_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    n = int(d["n"])
    ref = sym_onepol_C_shor(n, SUPPORT_CR0014, branch=POL_BRANCH)
    C = float(ref["C_shor_sym"])
    Cd = float(json.loads(Path("conjectures/active/L-0067.json").read_text())["C_dagger"]) if Path("conjectures/active/L-0067.json").is_file() else 9.562
    checks = {
        "support": list(d["support_cr0014"]) == list(SUPPORT_CR0014),
        "extends_cr0013": set(SUPPORT_GREEDY_SECOND).issubset(set(SUPPORT_CR0014)),
        "C_le_Cdagger": C <= Cd + 1e-9,
        "closes": d["closes_cr0014"] is True,
    }
    return all(checks.values()), checks


def save_l0068_certificate(
    cert: CertificateL0068,
    path: str | Path = "certificates/CERT-L0068-onepol-greedy-second-r24-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
