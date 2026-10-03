"""N5 certificate for L-0053 one-pol SOS Gram on C-R-0008."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0053_onepol_sos_gram_n5 import (
    C_FULLSYM_ONEPOL,
    D_ONEPOL,
    GRAM_JSON,
    lemma_l0053,
)
from ns_exploration.conjectures.l0052_sos_gram_n5 import (
    certified_C_ub_hi,
    load_gram_export,
    verify_gram_blocks,
)


@dataclass
class L0053OnePolSosGramCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    subclass: str
    radii: list[int]
    D: int
    C_dagger: float
    C_fullsym: float
    C_ub_float: float
    C_ub_hi: float
    ub_raw_hi: float
    gram_denom: int
    n_gram_blocks: int
    min_psd_eig: float
    max_rationalization_error: float
    beats_fullsym: bool
    closes_vs_Cdagger: bool
    closes_all_ic: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0053_certificate(
    n: int = 24,
    gram_path: str | Path = GRAM_JSON,
) -> L0053OnePolSosGramCertificate:
    b = lemma_l0053(n=n, gram_path=Path(gram_path))
    gram = load_gram_export(Path(gram_path))
    payload = {
        "n": n,
        "D": b.D,
        "subclass": b.subclass,
        "radii": b.radii,
        "C_ub_hi": b.C_ub_hi,
        "gram_denom": b.gram_denom,
        "n_gram_blocks": b.n_gram_blocks,
        "min_psd_eig": b.min_psd_eig,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0053OnePolSosGramCertificate(
        certificate_id="CERT-L0053-sos-gram-onepol-cr0008-N24",
        lemma_ref="L-0053 TSSOS order-2 + rational Gram on one-pol C-R-0008",
        route="B",
        evidence_level="N5",
        n=n,
        subclass=b.subclass,
        radii=list(b.radii or []),
        D=b.D,
        C_dagger=b.C_dagger,
        C_fullsym=b.C_fullsym,
        C_ub_float=b.C_ub_float,
        C_ub_hi=b.C_ub_hi,
        ub_raw_hi=b.ub_raw_hi,
        gram_denom=b.gram_denom,
        n_gram_blocks=b.n_gram_blocks,
        min_psd_eig=b.min_psd_eig,
        max_rationalization_error=b.max_rationalization_error,
        beats_fullsym=b.beats_fullsym,
        closes_vs_Cdagger=b.closes_vs_Cdagger,
        closes_all_ic=False,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0053_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    gram = load_gram_export(GRAM_JSON)
    gram_checks = verify_gram_blocks(gram)
    ub_raw = float(gram["ub_raw_abs"])
    _, C_hi = certified_C_ub_hi(ub_raw)
    checks = {
        **gram_checks,
        "D_matches": int(d["D"]) == D_ONEPOL,
        "C_ub_hi_valid": float(d["C_ub_hi"]) >= C_hi - 1e-12,
        "closes_vs_Cdagger": float(d["C_ub_hi"]) < float(d["C_dagger"]) - 1e-9,
        "beats_fullsym": float(d["C_ub_hi"]) < float(d["C_fullsym"]) - 1e-9,
        "C_fullsym_matches": abs(float(d["C_fullsym"]) - C_FULLSYM_ONEPOL) < 1e-6,
        "closes_all_ic_false": d["closes_all_ic"] is False,
        "evidence_level_N5": d.get("evidence_level") == "N5",
        "no_clay_claim": "None" in d.get("clay_implication", ""),
    }
    return all(checks.values()), checks


def save_l0053_certificate(
    cert: L0053OnePolSosGramCertificate,
    path: str | Path = "certificates/CERT-L0053-sos-gram-onepol-cr0008-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
