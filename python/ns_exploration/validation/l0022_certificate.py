"""N5 certificate for L-0022 gamma-bootstrap arithmetic."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0022_gamma_bootstrap import gamma_crit, lemma_l0022


@dataclass
class L0022GammaCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    E0: float
    nu: float
    T: float
    c0005_M: float
    Ncrit_c0005: float
    gamma_crit: float
    N_emp_max: float
    gamma_emp_max: float
    Omega_T_emp_max: float
    uniform_hypothesis_holds_empirically: bool
    proved_closes_c0005: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0022_certificate(
    n: int = 24,
    E0: float = 0.5,
    nu: float = 0.1,
    T: float = 0.02,
    empirical: bool = True,
    n_random: int = 8,
    ascent_steps: int = 15,
) -> L0022GammaCertificate:
    b = lemma_l0022(
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
        "c0005_M": b.c0005_M,
        "Ncrit_c0005": b.Ncrit_c0005,
        "gamma_crit": b.gamma_crit,
        "proved_closes_c0005": False,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()
    return L0022GammaCertificate(
        certificate_id="CERT-L0022-gamma-bootstrap-N24",
        lemma_ref="L-0022 gamma/N_* bootstrap arithmetic for C-0005 via L-0020",
        route="B",
        evidence_level="N5",
        n=n,
        E0=E0,
        nu=nu,
        T=T,
        c0005_M=b.c0005_M,
        Ncrit_c0005=b.Ncrit_c0005,
        gamma_crit=b.gamma_crit,
        N_emp_max=b.N_emp_max,
        gamma_emp_max=b.gamma_emp_max,
        Omega_T_emp_max=b.Omega_T_emp_max,
        uniform_hypothesis_holds_empirically=b.uniform_hypothesis_holds_empirically,
        proved_closes_c0005=False,
        clay_implication=b.clay_implication,
        payload_sha256=digest,
        notes=b.notes,
    )


def verify_l0022_certificate(d: dict) -> tuple[bool, dict[str, bool]]:
    from ns_exploration.conjectures.l0020_duhamel_h1 import lemma_l0020

    b20 = lemma_l0020(
        n=int(d["n"]),
        E0=float(d["E0"]),
        nu=float(d["nu"]),
        T=float(d["T"]),
        c0005_M=float(d["c0005_M"]),
    )
    g2 = gamma_crit(b20.E0, b20.c0005_M, b20.Ncrit_c0005)
    g = gamma_crit(float(d["E0"]), float(d["c0005_M"]), float(d["Ncrit_c0005"]))
    checks = {
        "Ncrit": abs(float(d["Ncrit_c0005"]) - b20.Ncrit_c0005) < 1e-8,
        "gamma_crit": abs(float(d["gamma_crit"]) - g2) < 1e-8,
        "gamma_matches_formula": abs(float(d["gamma_crit"]) - g) < 1e-10,
        "not_proved_closed": d["proved_closes_c0005"] is False,
    }
    return all(checks.values()), checks


def save_l0022_certificate(
    cert: L0022GammaCertificate,
    path: str | Path = "certificates/CERT-L0022-gamma-bootstrap-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
