"""N5 certificate for L-0023 cubic-dissipation C_★ arithmetic."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0023_cubic_dissipation import (
    C_star_for_M,
    lemma_l0023,
    omega_ode_bound,
)


@dataclass
class L0023CubicCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    E0: float
    nu: float
    T: float
    K2: int
    Omega0_cap: float
    c0005_M: float
    C_star: float
    Omega_at_Cstar: float
    C_emp_max: float
    emp_closes: bool
    proved_closes_c0005: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0023_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    empirical: bool = True,
    n_random: int = 16,
    ascent_steps: int = 12,
) -> L0023CubicCertificate:
    b = lemma_l0023(
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        empirical=empirical,
        n_random=n_random,
        ascent_steps=ascent_steps,
    )
    payload = {
        "n": n,
        "E0": E0,
        "nu": nu,
        "T": T,
        "K2": b.K2,
        "Omega0_cap": b.Omega0_cap,
        "c0005_M": b.c0005_M,
        "C_star": b.C_star,
        "Omega_at_Cstar": b.Omega_at_Cstar,
        "proved_closes_c0005": False,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0023CubicCertificate(
        certificate_id="CERT-L0023-cubic-dissipation-N24",
        lemma_ref="L-0023 cubic-dissipation ODE C_★ for C-0005",
        route="B",
        evidence_level="N5",
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        K2=b.K2,
        Omega0_cap=b.Omega0_cap,
        c0005_M=b.c0005_M,
        C_star=b.C_star,
        Omega_at_Cstar=b.Omega_at_Cstar,
        C_emp_max=b.C_emp_max,
        emp_closes=b.emp_closes,
        proved_closes_c0005=False,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0023_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    E0 = float(d["E0"])
    nu = float(d["nu"])
    T = float(d["T"])
    Omega0 = float(d["Omega0_cap"])
    M = float(d["c0005_M"])
    Cstar = C_star_for_M(E0, nu, T, Omega0, M)
    om = omega_ode_bound(float(d["C_star"]), E0, nu, T, Omega0)
    om_star = omega_ode_bound(Cstar, E0, nu, T, Omega0)
    checks = {
        "C_star": abs(float(d["C_star"]) - Cstar) < 1e-6,
        "Omega_at_Cstar_matches": abs(float(d["Omega_at_Cstar"]) - om_star) < 1e-4,
        "bound_at_cert_C_le_M": om <= M + 1e-6,
        "not_proved_closed": d["proved_closes_c0005"] is False,
        "K2_Omega0": abs(Omega0 - float(d["K2"]) * E0) < 1e-12,
    }
    return all(checks.values()), checks


def save_l0023_certificate(
    cert: L0023CubicCertificate,
    path: str | Path = "certificates/CERT-L0023-cubic-dissipation-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
