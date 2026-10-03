"""N5 certificate for L-0024 spectral-defect majorant / C-0005."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0024_spectral_defect import (
    dealias_shell_stats,
    lemma_l0024,
    max_shell_has_no_self_triads,
    worst_omega_T,
)


@dataclass
class L0024DefectCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    E0: float
    nu: float
    T: float
    K2: int
    m_K: int
    gap: int
    r_next: int
    rho_star: float
    max_shell_self_triads: int
    c0005_M: float
    Omega_T_worst: float
    Omega0_worst: float
    closes_c0005: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0024_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
) -> L0024DefectCertificate:
    b = lemma_l0024(n=n, E0=E0, nu=nu, T=T)
    payload = {
        "n": n,
        "E0": E0,
        "nu": nu,
        "T": T,
        "K2": b.K2,
        "m_K": b.m_K,
        "gap": b.gap,
        "rho_star": b.rho_star,
        "Omega_T_worst": b.Omega_T_worst,
        "closes_c0005": b.closes_c0005,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0024DefectCertificate(
        certificate_id="CERT-L0024-spectral-defect-C0005-N24",
        lemma_ref="L-0024 spectral-defect ODE closes C-0005",
        route="B",
        evidence_level="N5",
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=b.K2,
        m_K=b.m_K,
        gap=b.gap,
        r_next=b.r_next,
        rho_star=b.rho_star,
        max_shell_self_triads=0,
        c0005_M=b.c0005_M,
        Omega_T_worst=b.Omega_T_worst,
        Omega0_worst=b.Omega0_worst,
        closes_c0005=b.closes_c0005,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0024_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    stats = dealias_shell_stats(int(d["n"]))
    no_tri = max_shell_has_no_self_triads(int(d["n"]))
    omT, om0 = worst_omega_T(
        float(d["E0"]), float(d["nu"]), float(d["T"]), stats, n_grid=61
    )
    checks = {
        "no_self_triads": no_tri and int(d["max_shell_self_triads"]) == 0,
        "K2": int(d["K2"]) == stats["K2"],
        "m_K": int(d["m_K"]) == stats["m_K"],
        "gap": int(d["gap"]) == stats["gap"],
        "Omega_T": abs(float(d["Omega_T_worst"]) - omT) < 1e-4,
        "closes": bool(d["closes_c0005"]) and omT <= float(d["c0005_M"]) + 1e-8,
    }
    return all(checks.values()), checks


def save_l0024_certificate(
    cert: L0024DefectCertificate,
    path: str | Path = "certificates/CERT-L0024-spectral-defect-C0005-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
