"""N5 certificate for L-0019 Stokes floor constants."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0019_stokes_majorant import stokes_S
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class L0019StokesCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    E0: float
    nu: float
    t: float
    K2_star: int
    S_t: float
    Stokes_floor: float
    envelope_cap: float
    c0005_M: float
    floor_lt_c0005: bool
    floor_lt_envelope: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0019_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    t: float = 0.02,
    c0005_M: float = 61.23693461895651,
) -> L0019StokesCertificate:
    S, r_star = stokes_S(n, nu, t)
    floor = S * E0
    K2, _ = full_dealias_exact_stats(n)
    env = float(K2) * E0
    payload = {
        "n": n,
        "E0": E0,
        "nu": nu,
        "t": t,
        "K2_star": r_star,
        "S_t": S,
        "Stokes_floor": floor,
        "envelope_cap": env,
        "c0005_M": c0005_M,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0019StokesCertificate(
        certificate_id="CERT-L0019-stokes-floor-N24",
        lemma_ref="L-0019 Stokes spectral majorant (linear)",
        route="B",
        evidence_level="N5",
        n=n,
        E0=E0,
        nu=nu,
        t=t,
        K2_star=r_star,
        S_t=S,
        Stokes_floor=floor,
        envelope_cap=env,
        c0005_M=c0005_M,
        floor_lt_c0005=floor < c0005_M,
        floor_lt_envelope=floor < env,
        clay_implication="None. Linear Stokes only; not nonlinear NS; not Clay.",
        payload_sha256=digest,
        notes=f"K²_star={r_star}, floor={floor:.12g}, gap_to_C0005={c0005_M - floor:.6g}.",
    )


def verify_l0019_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    S, r = stokes_S(int(d["n"]), float(d["nu"]), float(d["t"]))
    floor = S * float(d["E0"])
    K2, _ = full_dealias_exact_stats(int(d["n"]))
    env = float(K2) * float(d["E0"])
    checks = {
        "K2_star": int(d["K2_star"]) == r,
        "S_match": abs(float(d["S_t"]) - S) < 1e-12,
        "floor_match": abs(float(d["Stokes_floor"]) - floor) < 1e-12,
        "envelope_match": abs(float(d["envelope_cap"]) - env) < 1e-12,
        "floor_lt_c0005": float(d["Stokes_floor"]) < float(d["c0005_M"]),
        "floor_lt_envelope": float(d["Stokes_floor"]) < float(d["envelope_cap"]),
    }
    return all(checks.values()), checks


def save_l0019_certificate(
    cert: L0019StokesCertificate,
    path: str | Path = "certificates/CERT-L0019-stokes-floor-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
