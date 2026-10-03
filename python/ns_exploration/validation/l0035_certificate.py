"""Certificate for L-0035 shell two-point ‖∇u‖_∞ stretch techo."""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0035_ginf_stretch import (
    C_eff_ginf,
    lemma_l0035,
    shell_radii_masses,
)


@dataclass
class CertificateL0035:
    cert_id: str
    lemma_id: str
    n: int
    c0007_M: float
    C_dagger: float
    C_at_equipartition: float
    C_max_low_slab: float
    Omega_c: float
    ginf_ode_worst: float
    minprod_ode_worst: float
    closes_c0007: bool
    clay_implication: str

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0035_certificate(n: int = 24) -> CertificateL0035:
    b = lemma_l0035(n=n, n_grid=7)
    return CertificateL0035(
        cert_id="CERT-L0035-ginf-stretch-C0007-N24",
        lemma_id="L-0035",
        n=b.n,
        c0007_M=b.c0007_M,
        C_dagger=b.C_dagger,
        C_at_equipartition=b.C_at_equipartition,
        C_max_low_slab=b.C_max_low_slab,
        Omega_c=b.Omega_c,
        ginf_ode_worst=b.ginf_ode_worst,
        minprod_ode_worst=b.minprod_ode_worst,
        closes_c0007=b.closes_c0007,
        clay_implication=b.clay_implication,
    )


def verify_l0035_certificate(d: dict) -> tuple[bool, dict]:
    radii, masses = shell_radii_masses(int(d["n"]))
    C_eq = C_eff_ginf(0.5, 0.5, radii, masses)
    checks = {
        "C_eq_matches": abs(C_eq - float(d["C_at_equipartition"])) < 1e-6,
        "C_eq_le_Cdagger": float(d["C_at_equipartition"]) <= float(d["C_dagger"]),
        "C_max_above_Cdagger": float(d["C_max_low_slab"]) > float(d["C_dagger"]),
        "Omega_c_near_E0": abs(float(d["Omega_c"]) - 0.5) < 0.05,
        "ginf_ode_above_M": float(d["ginf_ode_worst"]) > float(d["c0007_M"]),
        "minprod_ode_above_M": float(d["minprod_ode_worst"]) > float(d["c0007_M"]),
        "closes_false": d["closes_c0007"] is False,
    }
    return all(checks.values()), checks


def save_l0035_certificate(
    cert: CertificateL0035,
    path: str | Path = "certificates/CERT-L0035-ginf-stretch-C0007-N24.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
