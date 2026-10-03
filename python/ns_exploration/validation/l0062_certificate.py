"""Certificates for L-0062 COSMO band ladder."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0052_sos_gram_n5 import certified_C_ub_hi
from ns_exploration.conjectures.l0062_cosmo_ladder_n5 import COSMO_BANDS, lemma_l0062


@dataclass
class L0062BandCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    band_radii: list[int]
    D: int
    C_dagger: float
    C_fullsym: float
    C_ub_float: float
    C_ub_hi: float
    ub_raw_hi: float
    gram_denom: int
    n_gram_blocks: int
    beats_fullsym: bool
    closes_vs_Cdagger: bool
    closes_all_ic: bool
    solver: str
    related_cr: str
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0062_certificates(n: int = 24) -> list[L0062BandCertificate]:
    bound = lemma_l0062(n=n)
    certs: list[L0062BandCertificate] = []
    for row, spec in zip(bound.rows or [], COSMO_BANDS, strict=True):
        payload = {
            "n": n,
            "D": row.D,
            "radii": row.band_radii,
            "C_ub_hi": row.C_ub_hi,
        }
        digest = hashlib.sha256(
            json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        certs.append(
            L0062BandCertificate(
                certificate_id=row.cert_id,
                lemma_ref=f"L-0062 N5 interval COSMO Hermitian {{1..{row.band_radii[-1]}}}",
                route="B",
                evidence_level="N5",
                n=n,
                band_radii=row.band_radii,
                D=row.D,
                C_dagger=row.C_dagger,
                C_fullsym=row.C_fullsym,
                C_ub_float=row.C_ub_float,
                C_ub_hi=row.C_ub_hi,
                ub_raw_hi=row.ub_raw_hi,
                gram_denom=0,
                n_gram_blocks=row.gram_blocks,
                beats_fullsym=row.C_ub_hi < row.C_fullsym - 1e-9,
                closes_vs_Cdagger=row.C_ub_hi < row.C_dagger - 1e-9,
                closes_all_ic=False,
                solver="COSMO",
                related_cr=row.related_cr,
                clay_implication="None. Band subclass only; not all-IC; not Clay.",
                payload_sha256=digest,
                notes=row.notes,
            )
        )
    return certs


def verify_l0062_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    spec = next(s for s in COSMO_BANDS if s["cert_id"] == d["certificate_id"])
    tssos = json.loads((Path(spec["data_dir"]) / spec["tssos_file"]).read_text(encoding="utf-8"))
    ub_raw = float(tssos["ub_raw_abs"])
    _, C_hi = certified_C_ub_hi(ub_raw)
    checks = {
        "D_matches": int(d["D"]) == int(spec["D"]),
        "C_ub_hi_valid": float(d["C_ub_hi"]) >= C_hi - 1e-12,
        "closes_vs_Cdagger": float(d["C_ub_hi"]) < float(d["C_dagger"]) - 1e-9,
        "beats_fullsym": float(d["C_ub_hi"]) < float(d["C_fullsym"]) - 1e-9,
        "closes_all_ic_false": d["closes_all_ic"] is False,
        "evidence_level_N5": d.get("evidence_level") == "N5",
    }
    return all(checks.values()), checks


def save_l0062_certificates(
    certs: list[L0062BandCertificate],
    cert_dir: str | Path = "certificates",
) -> list[Path]:
    cert_dir = Path(cert_dir)
    cert_dir.mkdir(parents=True, exist_ok=True)
    paths: list[Path] = []
    for cert in certs:
        p = cert_dir / f"{cert.certificate_id}.json"
        p.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
        paths.append(p)
    return paths
