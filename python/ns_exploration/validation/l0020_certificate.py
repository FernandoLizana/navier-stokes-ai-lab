"""N5 certificate for L-0020 Stokes–Duhamel H¹ constants."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0020_duhamel_h1 import (
    Ncrit_for_M,
    full_mask_rho_star,
    integrate_sigma,
    lemma_l0020,
    omega_duhamel_H1,
    young_Nstar,
)
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class L0020DuhamelCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    E0: float
    nu: float
    T: float
    K2: int
    rho_star: float
    S_T: float
    I_sigma: float
    Nstar_young: float
    Omega_young_H1: float
    envelope_cap: float
    c0005_M: float
    Ncrit_c0005: float
    young_closes_c0005: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0020_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    c0005_M: float = 61.23693461895651,
) -> L0020DuhamelCertificate:
    b = lemma_l0020(n=n, E0=E0, nu=nu, T=T, c0005_M=c0005_M)
    payload = {
        "n": n,
        "E0": E0,
        "nu": nu,
        "T": T,
        "K2": b.K2,
        "rho_star": b.rho_star,
        "S_T": b.S_T,
        "I_sigma": b.I_sigma,
        "Nstar_young": b.Nstar_young,
        "Omega_young_H1": b.Omega_young_H1,
        "Ncrit_c0005": b.Ncrit_c0005,
        "c0005_M": c0005_M,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0020DuhamelCertificate(
        certificate_id="CERT-L0020-duhamel-H1-N24",
        lemma_ref="L-0020 Stokes–Duhamel H¹ majorant",
        route="B",
        evidence_level="N5",
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=b.K2,
        rho_star=b.rho_star,
        S_T=b.S_T,
        I_sigma=b.I_sigma,
        Nstar_young=b.Nstar_young,
        Omega_young_H1=b.Omega_young_H1,
        envelope_cap=b.envelope_cap,
        c0005_M=c0005_M,
        Ncrit_c0005=b.Ncrit_c0005,
        young_closes_c0005=b.young_closes_c0005,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0020_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    n = int(d["n"])
    E0, nu, T = float(d["E0"]), float(d["nu"]), float(d["T"])
    cM = float(d["c0005_M"])
    b = lemma_l0020(n=n, E0=E0, nu=nu, T=T, c0005_M=cM)
    K2, _ = full_dealias_exact_stats(n)
    rho, _, _ = full_mask_rho_star(n)
    checks = {
        "K2": int(d["K2"]) == int(K2),
        "rho": abs(float(d["rho_star"]) - rho) < 1e-9,
        "S_T": abs(float(d["S_T"]) - b.S_T) < 1e-10,
        "I_sigma": abs(float(d["I_sigma"]) - b.I_sigma) < 1e-8,
        "Nstar": abs(float(d["Nstar_young"]) - young_Nstar(rho, E0)) < 1e-10,
        "Omega_formula": abs(
            float(d["Omega_young_H1"])
            - omega_duhamel_H1(b.S_T, E0, b.Nstar_young, b.I_sigma)
        )
        < 1e-8,
        "Ncrit": abs(float(d["Ncrit_c0005"]) - Ncrit_for_M(b.S_T, E0, b.I_sigma, cM))
        < 1e-8,
        "young_flag": bool(d["young_closes_c0005"]) == (b.Omega_young_H1 <= cM + 1e-12),
    }
    return all(checks.values()), checks


def save_l0020_certificate(
    cert: L0020DuhamelCertificate,
    path: str | Path = "certificates/CERT-L0020-duhamel-H1-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
