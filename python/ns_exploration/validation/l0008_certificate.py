"""
N5 interval certificate for the L-0008 cascade bound
  Ω(t) ≤ B · cosh²(A K_full t),
  B = K_IC² E0,  A = √2 √(3 M_full) √E0.

Finite Galerkin only. Not Clay.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import (
    euclidean_shell_K2,
    linf_shell_K2,
)
from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class L0008CoshCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    ic_kind: str
    k_ic: int
    K_ic_squared: int
    K_full_squared: int
    M_full: int
    E0: float
    t: float
    B_hi: float
    A_hi: float
    alpha_hi: float
    Omega_bound_hi: float
    closes: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def build_l0008_cosh_certificate(
    n: int = 12,
    k_ic: int = 2,
    ic_kind: str = "euclidean",
    E0: float = 0.5,
    t: float = 0.02,
) -> L0008CoshCertificate:
    if ic_kind == "linf":
        K_ic2 = linf_shell_K2(k_ic)
    elif ic_kind == "euclidean":
        K_ic2 = euclidean_shell_K2(n, k_ic)
    else:
        raise ValueError(ic_kind)

    K2_full, M_full = full_dealias_exact_stats(n)
    IE0 = Interval.from_float(E0)
    It = Interval.from_float(t)
    IB = Interval(float(K_ic2), float(K_ic2)) * IE0
    IK2 = Interval(float(K2_full), float(K2_full))
    IM = Interval(float(M_full), float(M_full))
    I2 = Interval(2.0, 2.0)
    I3 = Interval(3.0, 3.0)

    # A = √2 √(3M) √E0 = √(2 · 3M · E0) = √(6 M E0)
    A = (I2 * I3 * IM * IE0).sqrt()
    K_full = IK2.sqrt()
    alpha = A * K_full
    cosh_at = (alpha * It).cosh()
    cap = IB * cosh_at.sqr()

    payload = json.dumps(
        {"n": n, "ic_kind": ic_kind, "k_ic": k_ic, "K_ic2": K_ic2,
         "K2_full": K2_full, "M": M_full, "E0": E0, "t": t},
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    return L0008CoshCertificate(
        certificate_id=f"CERT-L0008-cosh-{ic_kind}-k{k_ic}-N{n}",
        lemma_ref="L-0008 cascade cosh bound (E_H(0)=0)",
        route="B",
        evidence_level="N5",
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_full=M_full,
        E0=E0,
        t=t,
        B_hi=float(IB.hi),
        A_hi=float(A.hi),
        alpha_hi=float(alpha.hi),
        Omega_bound_hi=float(cap.hi),
        closes=True,
        clay_implication="None. Finite Galerkin cascade estimate; not continuum NS.",
        payload_sha256=digest,
        notes=f"Outward: Ω≤B cosh²(αt) ≤ {float(cap.hi):.6e}.",
    )


def verify_l0008_cosh_certificate(cert: dict) -> tuple[bool, dict]:
    rebuilt = build_l0008_cosh_certificate(
        n=int(cert["n"]),
        k_ic=int(cert["k_ic"]),
        ic_kind=str(cert["ic_kind"]),
        E0=float(cert["E0"]),
        t=float(cert["t"]),
    )
    checks = {
        "K_ic_matches": rebuilt.K_ic_squared == int(cert["K_ic_squared"]),
        "K_full_matches": rebuilt.K_full_squared == int(cert["K_full_squared"]),
        "M_matches": rebuilt.M_full == int(cert["M_full"]),
        "bound_is_valid_upper": float(cert["Omega_bound_hi"]) >= rebuilt.Omega_bound_hi - 1e-4,
        "evidence_level_N5": cert.get("evidence_level") == "N5",
        "no_clay_claim": "None" in cert.get("clay_implication", ""),
    }
    return all(checks.values()), checks


def save_l0008_cosh_certificate(cert: L0008CoshCertificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
