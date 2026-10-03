"""N5 certificate for L-0021 mono-radial quartic Gram + Duhamel."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0021_quartic_shell import lemma_l0021


@dataclass
class L0021QuarticCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    E0: float
    nu: float
    T: float
    alpha_star: float
    r_star: int
    n_shells: int
    Omega_radial_H1: float
    envelope_cap: float
    beats_envelope: bool
    radial_M: float
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0021_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    iters: int = 50,
) -> L0021QuarticCertificate:
    b = lemma_l0021(n=n, E0=E0, nu=nu, T=T, iters=iters)
    radial_M = float(b.Omega_radial_H1)
    payload = {
        "n": n,
        "E0": E0,
        "nu": nu,
        "T": T,
        "alpha_star": b.alpha_star,
        "r_star": b.r_star,
        "n_shells": b.n_shells,
        "Omega_radial_H1": b.Omega_radial_H1,
        "envelope_cap": b.envelope_cap,
        "radial_M": radial_M,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0021QuarticCertificate(
        certificate_id="CERT-L0021-quartic-radial-N24",
        lemma_ref="L-0021 mono-radial quartic Gram + L-0020 Duhamel",
        route="B",
        evidence_level="N5",
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        alpha_star=b.alpha_star,
        r_star=b.r_star,
        n_shells=b.n_shells,
        Omega_radial_H1=b.Omega_radial_H1,
        envelope_cap=b.envelope_cap,
        beats_envelope=b.beats_envelope,
        radial_M=radial_M,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0021_certificate(d: dict, iters: int = 50) -> tuple[bool, dict[str, bool]]:
    b = lemma_l0021(
        n=int(d["n"]),
        E0=float(d["E0"]),
        nu=float(d["nu"]),
        T=float(d["T"]),
        iters=iters,
    )
    checks = {
        "alpha_star": abs(float(d["alpha_star"]) - b.alpha_star) < 1e-6,
        "r_star": int(d["r_star"]) == b.r_star,
        "n_shells": int(d["n_shells"]) == b.n_shells,
        "Omega": abs(float(d["Omega_radial_H1"]) - b.Omega_radial_H1) < 1e-6,
        "beats_envelope": bool(d["beats_envelope"]) == (b.Omega_radial_H1 < b.envelope_cap),
        "radial_M": abs(float(d["radial_M"]) - b.Omega_radial_H1) < 1e-9,
    }
    return all(checks.values()), checks


def save_l0021_certificate(
    cert: L0021QuarticCertificate,
    path: str | Path = "certificates/CERT-L0021-quartic-radial-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
