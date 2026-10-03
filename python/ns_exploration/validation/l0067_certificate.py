"""N5 certificate for L-0067 one-pol greedy-second + shell 24."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0052_sos_gram_n5 import certified_C_ub_hi
from ns_exploration.conjectures.l0067_onepol_greedy_r24_sos import SUPPORT_CR0014, lemma_l0067

CERT_ID = "CERT-L0067-sos-onepol-greedy-second-r24-N24"


@dataclass
class L0067Certificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    support_radii: list[int]
    pol_branch: str
    D: int
    C_dagger: float
    C_shor_sym: float
    C_ub_float: float
    C_ub_hi: float
    ub_raw_hi: float
    closes_vs_Cdagger: bool
    proposed_cr: str
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0067_certificate(n: int = 24) -> L0067Certificate:
    b = lemma_l0067(n=n)
    payload = {"n": n, "D": b.D, "support": b.support_radii, "C_ub_hi": b.C_ub_hi}
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0067Certificate(
        certificate_id=CERT_ID,
        lemma_ref="L-0067 N5 interval COSMO one-pol 12 shells",
        route="B",
        evidence_level="N5",
        n=n,
        support_radii=list(SUPPORT_CR0014),
        pol_branch=b.pol_branch,
        D=b.D,
        C_dagger=b.C_dagger,
        C_shor_sym=b.C_shor_sym,
        C_ub_float=b.C_ub_float,
        C_ub_hi=b.C_ub_hi,
        ub_raw_hi=b.ub_raw_hi,
        closes_vs_Cdagger=b.closes_vs_Cdagger,
        proposed_cr="C-R-0014",
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0067_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    tssos = json.loads(
        Path("tools/sos_julia/data/onepol_greedy_second_r24/cosmo_tssos_result.json").read_text(
            encoding="utf-8"
        )
    )
    _, C_hi = certified_C_ub_hi(float(tssos["ub_raw_abs"]))
    checks = {
        "C_ub_hi_matches": abs(C_hi - float(d["C_ub_hi"])) < 1e-12,
        "closes_vs_Cdagger": float(d["C_ub_hi"]) < float(d["C_dagger"]) - 1e-9,
    }
    return all(checks.values()), checks


def save_l0067_certificate(
    cert: L0067Certificate,
    path: str | Path = f"certificates/{CERT_ID}.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
