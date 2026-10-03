"""
N5 certificate for the L-0010 short-time C-S bound:

  On N=16, IC |k|_∞≤4, E0=1/2, for all t ≤ T_cert,
  Ω(t) ≤ R  with R = C-S-0002 proposed M.

Verified by outward interval check of the Duhamel bootstrap inequality
  T² (U_L + W √ε)² R ≤ ε,   ε = (R-B)/K_full².
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2
from ns_exploration.conjectures.l0009_twoscale import shell_mode_count
from ns_exploration.conjectures.l0010_duhamel_short import max_T_duhamel
from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class L0010ShortCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    ic_kind: str
    k_ic: int
    K_ic_squared: int
    K_full_squared: int
    M_L: int
    M_H: int
    E0: float
    R: float
    T_cert: float
    T_star: float
    B_hi: float
    eps_lo: float
    E_H_cap_hi: float
    bootstrap_holds: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _interval_bootstrap(
    K_ic2: int,
    K2_full: int,
    M_L: int,
    M_H: int,
    E0: float,
    R: float,
    T: float,
) -> tuple[Interval, Interval, bool]:
    IE0 = Interval.from_float(E0)
    IR = Interval.from_float(R)
    IT = Interval.from_float(T)
    IB = Interval(float(K_ic2), float(K_ic2)) * IE0
    IK2 = Interval(float(K2_full), float(K2_full))
    IUL = (Interval(3.0, 3.0) * Interval(float(M_L), float(M_L))).sqrt() * (
        Interval(2.0, 2.0) * IE0
    ).sqrt()
    IW = (
        (Interval(3.0, 3.0) * Interval(float(M_H), float(M_H))).sqrt()
        * Interval(2.0, 2.0).sqrt()
        if M_H > 0
        else Interval(0.0, 0.0)
    )
    diff = IR - IB
    if diff.lo <= 0:
        raise RuntimeError("R must exceed B")
    eps = diff / IK2
    E_H_cap = IT.sqr() * (IUL + IW * eps.sqrt()).sqr() * IR
    holds = bool(E_H_cap.hi <= eps.lo + 1e-12)
    return eps, E_H_cap, holds


def build_l0010_short_certificate(
    n: int = 16,
    k_ic: int = 4,
    ic_kind: str = "linf",
    E0: float = 0.5,
    R: float = 48.15928260103266,
    safety: float = 0.99,
) -> L0010ShortCertificate:
    if ic_kind != "linf":
        raise ValueError("certificate implemented for linf C-S class")

    K_ic2 = linf_shell_K2(k_ic)
    K2_full, M_full = full_dealias_exact_stats(n)
    M_L = shell_mode_count(n, k_ic, ic_kind)
    M_H = M_full - M_L
    K_full = math.sqrt(float(K2_full))
    B = float(K_ic2) * E0
    U_L = math.sqrt(3.0 * M_L) * math.sqrt(2.0 * E0)
    W = math.sqrt(3.0 * M_H) * math.sqrt(2.0) if M_H > 0 else 0.0

    T_star = max_T_duhamel(B, K_full, U_L, W, R)
    T_cert = float(safety) * T_star
    eps, E_H_cap, holds = _interval_bootstrap(K_ic2, K2_full, M_L, M_H, E0, R, T_cert)
    if not holds:
        raise RuntimeError("bootstrap failed at T_cert; lower safety")

    payload = json.dumps(
        {"n": n, "k_ic": k_ic, "R": R, "T_cert": T_cert, "M_L": M_L, "M_H": M_H},
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    return L0010ShortCertificate(
        certificate_id="CERT-L0010-duhamel-CS0002-short-N16",
        lemma_ref="L-0010 Duhamel bootstrap (short-time C-S-0002 bound)",
        route="B",
        evidence_level="N5",
        n=n,
        ic_kind=ic_kind,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        E0=E0,
        R=R,
        T_cert=T_cert,
        T_star=T_star,
        B_hi=float((Interval(float(K_ic2), float(K_ic2)) * Interval.from_float(E0)).hi),
        eps_lo=float(eps.lo),
        E_H_cap_hi=float(E_H_cap.hi),
        bootstrap_holds=holds,
        clay_implication=(
            "None. Short-time finite-Galerkin bound only; does not prove C-S-0002 at T=0.02."
        ),
        payload_sha256=digest,
        notes=(
            f"T_cert={T_cert:.8f} (= {safety} * T_star={T_star:.8f}). "
            f"E_H_cap_hi={float(E_H_cap.hi):.8e} <= eps_lo={float(eps.lo):.8e}."
        ),
    )


def verify_l0010_short_certificate(cert: dict) -> tuple[bool, dict]:
    K_ic2 = linf_shell_K2(int(cert["k_ic"]))
    K2_full, M_full = full_dealias_exact_stats(int(cert["n"]))
    M_L = shell_mode_count(int(cert["n"]), int(cert["k_ic"]), "linf")
    M_H_re = M_full - M_L

    eps, E_H_cap, holds = _interval_bootstrap(
        K_ic2,
        K2_full,
        M_L,
        M_H_re,
        float(cert["E0"]),
        float(cert["R"]),
        float(cert["T_cert"]),
    )
    checks = {
        "K_ic_matches": K_ic2 == int(cert["K_ic_squared"]),
        "K_full_matches": K2_full == int(cert["K_full_squared"]),
        "M_L_matches": M_L == int(cert["M_L"]),
        "M_H_matches": M_H_re == int(cert["M_H"]),
        "bootstrap_holds": holds,
        "stored_flag": bool(cert.get("bootstrap_holds")),
        "E_H_cap_valid_upper": float(cert["E_H_cap_hi"]) >= float(E_H_cap.hi) - 1e-9,
        "evidence_level_N5": cert.get("evidence_level") == "N5",
        "no_clay_claim": "None" in cert.get("clay_implication", ""),
        "short_horizon_only": float(cert["T_cert"]) < 0.02,
    }
    return all(checks.values()), checks


def save_l0010_short_certificate(cert: L0010ShortCertificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
