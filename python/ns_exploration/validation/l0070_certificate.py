"""N5 certificate for L-0070 all-IC C-0007 closure via L-0024."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0070_c0007_all_ic_l0024 import lemma_l0070
from ns_exploration.conjectures.l0024_spectral_defect import (
    dealias_shell_stats,
    max_shell_has_no_self_triads,
    worst_omega_T,
)
from ns_exploration.conjectures.l0026_c0007_techo import l0024_monotone_in_omega0


@dataclass
class L0070AllIcCertificate:
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
    L0024_majorant: float
    proved_bound_M: float
    sharp_target_M: float
    Omega0_worst: float
    monotone_in_Omega0: bool
    closes_c0007_all_ic: bool
    closes_c0008_sharp: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0070_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
) -> L0070AllIcCertificate:
    b = lemma_l0070(n=n, E0=E0, nu=nu, T=T)
    payload = {
        "n": n,
        "E0": E0,
        "nu": nu,
        "T": T,
        "K2": b.K2,
        "proved_bound_M": b.proved_bound_M,
        "L0024_majorant": b.L0024_majorant,
        "closes_c0007_all_ic": b.closes_c0007_all_ic,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0070AllIcCertificate(
        certificate_id="CERT-L0070-all-ic-C0007-N24",
        lemma_ref="L-0070 all-IC C-0007 via L-0024 spectral-defect",
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
        L0024_majorant=b.L0024_majorant,
        proved_bound_M=b.proved_bound_M,
        sharp_target_M=b.sharp_target_M,
        Omega0_worst=b.Omega0_worst,
        monotone_in_Omega0=b.monotone_in_Omega0,
        closes_c0007_all_ic=b.closes_c0007_all_ic,
        closes_c0008_sharp=b.closes_c0008_sharp,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0070_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    stats = dealias_shell_stats(int(d["n"]))
    no_tri = max_shell_has_no_self_triads(int(d["n"]))
    omT, _ = worst_omega_T(
        float(d["E0"]), float(d["nu"]), float(d["T"]), stats, n_grid=61
    )
    monotone = l0024_monotone_in_omega0(
        float(d["E0"]), float(d["nu"]), float(d["T"]), stats, n_grid=61
    )
    M = float(d["proved_bound_M"])
    checks = {
        "no_self_triads": no_tri,
        "K2": int(d["K2"]) == stats["K2"],
        "monotone": monotone and bool(d["monotone_in_Omega0"]),
        "Omega_T": abs(float(d["L0024_majorant"]) - omT) < 1e-4,
        "closes_all_ic": bool(d["closes_c0007_all_ic"]) and omT <= M + 1e-8,
        "sharp_still_open": not bool(d["closes_c0008_sharp"]),
    }
    return all(checks.values()), checks


def save_l0070_certificate(
    cert: L0070AllIcCertificate,
    path: str | Path = "certificates/CERT-L0070-all-ic-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
