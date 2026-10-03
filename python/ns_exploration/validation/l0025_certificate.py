"""N5 certificate for L-0025 / C-0006 (N≤32 spectral defect)."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0024_spectral_defect import max_shell_has_no_self_triads
from ns_exploration.conjectures.l0025_multi_n_defect import lemma_l0025


@dataclass
class L0025MultiNCertificate:
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
    Stokes_floor: float
    c0006_M: float
    Omega_T_worst: float
    closes_c0006: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0025_certificate(
    n: int = 32,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
) -> L0025MultiNCertificate:
    b = lemma_l0025(n=n, E0=E0, nu=nu, T=T)
    payload = {
        "n": n,
        "E0": E0,
        "nu": nu,
        "T": T,
        "K2": b.K2,
        "gap": b.gap,
        "c0006_M": b.c0006_M,
        "Omega_T_worst": b.Omega_T_worst,
        "closes_c0006": b.closes_c0006,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0025MultiNCertificate(
        certificate_id="CERT-L0025-spectral-defect-C0006-N32",
        lemma_ref="L-0025 spectral-defect ODE closes C-0006 (N≤32)",
        route="B",
        evidence_level="N5",
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=b.K2,
        m_K=b.m_K,
        gap=b.gap,
        Stokes_floor=b.Stokes_floor,
        c0006_M=b.c0006_M,
        Omega_T_worst=b.Omega_T_worst,
        closes_c0006=b.closes_c0006,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0025_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    b = lemma_l0025(
        n=int(d["n"]),
        E0=float(d["E0"]),
        nu=float(d["nu"]),
        T=float(d["T"]),
        n_grid=41,
    )
    checks = {
        "no_self_triads": max_shell_has_no_self_triads(int(d["n"])),
        "K2": int(d["K2"]) == b.K2,
        "gap": int(d["gap"]) == b.gap,
        "Omega_T": abs(float(d["Omega_T_worst"]) - b.Omega_T_worst) < 1e-3,
        "closes": bool(d["closes_c0006"]) and b.closes_c0006,
        "M": abs(float(d["c0006_M"]) - b.c0006_M) < 1e-8,
    }
    return all(checks.values()), checks


def save_l0025_certificate(
    cert: L0025MultiNCertificate,
    path: str | Path = "certificates/CERT-L0025-spectral-defect-C0006-N32.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
