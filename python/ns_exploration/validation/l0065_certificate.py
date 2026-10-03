"""N5 certificate for L-0065 one-pol greedy-second SOS."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0052_sos_gram_n5 import certified_C_ub_hi
from ns_exploration.conjectures.l0065_onepol_greedy_second_sos import (
    SUPPORT_GREEDY_SECOND,
    lemma_l0065,
)

CERT_ID = "CERT-L0065-sos-onepol-greedy-second-N24"


@dataclass
class L0065Certificate:
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
    C_fullsym: float
    C_ub_float: float
    C_ub_hi: float
    ub_raw_hi: float
    beats_fullsym: bool
    beats_shor: bool
    closes_vs_Cdagger: bool
    closes_all_ic: bool
    solver: str
    proposed_cr: str
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0065_certificate(n: int = 24) -> L0065Certificate:
    b = lemma_l0065(n=n)
    payload = {
        "n": n,
        "D": b.D,
        "support": b.support_radii,
        "branch": b.pol_branch,
        "C_ub_hi": b.C_ub_hi,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0065Certificate(
        certificate_id=CERT_ID,
        lemma_ref="L-0065 N5 interval COSMO one-pol greedy-second",
        route="B",
        evidence_level="N5",
        n=n,
        support_radii=list(SUPPORT_GREEDY_SECOND),
        pol_branch=b.pol_branch,
        D=b.D,
        C_dagger=b.C_dagger,
        C_shor_sym=b.C_shor_sym,
        C_fullsym=b.C_fullsym,
        C_ub_float=b.C_ub_float,
        C_ub_hi=b.C_ub_hi,
        ub_raw_hi=b.ub_raw_hi,
        beats_fullsym=b.beats_fullsym,
        beats_shor=b.beats_shor,
        closes_vs_Cdagger=b.closes_vs_Cdagger,
        closes_all_ic=False,
        solver=b.solver,
        proposed_cr=b.proposed_cr,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0065_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    tssos = json.loads(
        Path("tools/sos_julia/data/onepol_greedy_second/cosmo_tssos_result.json").read_text(
            encoding="utf-8"
        )
    )
    ub_raw = float(tssos["ub_raw_abs"])
    _, C_hi = certified_C_ub_hi(ub_raw)
    checks = {
        "C_ub_hi_matches_interval": abs(C_hi - float(d["C_ub_hi"])) < 1e-12,
        "closes_vs_Cdagger": float(d["C_ub_hi"]) < float(d["C_dagger"]) - 1e-9,
        "beats_shor": float(d["C_ub_hi"]) < float(d["C_shor_sym"]) - 1e-9,
    }
    return all(checks.values()), checks


def save_l0065_certificate(
    cert: L0065Certificate,
    path: str | Path = f"certificates/{CERT_ID}.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
