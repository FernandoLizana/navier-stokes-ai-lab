"""
N5 certificate for L-0011 rigorous short-time C-S bound
(div-free embedding + viscous Duhamel bootstrap).
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2
from ns_exploration.conjectures.l0009_twoscale import shell_mode_count
from ns_exploration.conjectures.l0010_duhamel_short import kappa_H_min
from ns_exploration.conjectures.l0011_sharpened_short import (
    embedding_U,
    embedding_W,
    max_T_bootstrap,
    tau_viscous,
)
from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class L0011ShortCertificate:
    certificate_id: str
    lemma_ref: str
    route: str
    evidence_level: str
    n: int
    k_ic: int
    K_ic_squared: int
    K_full_squared: int
    M_L: int
    M_H: int
    kappa_H: float
    E0: float
    nu: float
    R: float
    T_cert: float
    T_star: float
    tau_hi: float
    eps_lo: float
    E_H_cap_hi: float
    bootstrap_holds: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _bootstrap_interval(
    K_ic2: int,
    K2_full: int,
    M_L: int,
    M_H: int,
    E0: float,
    nu: float,
    kappa_H: float,
    R: float,
    T: float,
) -> tuple[Interval, Interval, Interval, bool]:
    IE0 = Interval.from_float(E0)
    IR = Interval.from_float(R)
    IT = Interval.from_float(T)
    Inu = Interval.from_float(nu)
    Ikap = Interval.from_float(kappa_H)
    IB = Interval(float(K_ic2), float(K_ic2)) * IE0
    IK2 = Interval(float(K2_full), float(K2_full))
    # U_L = √(2 M_L) √(2 E0) = 2 √(M_L E0 / 1)? √(2M)√(2E)=2√(M E)
    IUL = (Interval(2.0, 2.0) * Interval(float(M_L), float(M_L))).sqrt() * (
        Interval(2.0, 2.0) * IE0
    ).sqrt()
    IW = (
        (Interval(2.0, 2.0) * Interval(2.0, 2.0) * Interval(float(M_H), float(M_H))).sqrt()
        if M_H > 0
        else Interval(0.0, 0.0)
    )  # √(2*2*M_H)=√(4 M_H)=2√M_H; embedding_W=√(2*2*M)=√(4M)=2√M. Good.

    # τ = (1 - exp(-ν κ T)) / (ν κ)
    a = Inu * Ikap
    # exp(-a T): use Interval.exp on (-a*T)
    neg_aT = Interval(0.0, 0.0) - (a * IT)
    # For upper bound of τ we need careful enclosure.
    # τ = (1-e^{-x})/ (x/T) wait x=aT, τ=(1-e^{-x})/a
    x = a * IT
    one = Interval(1.0, 1.0)
    # e^{-x}: if x≥0, e^{-x} ∈ [e^{-x.hi}, e^{-x.lo}]
    import numpy as np

    e_lo = float(np.nextafter(math.exp(-x.hi), -np.inf)) if x.hi >= 0 else float("nan")
    e_hi = float(np.nextafter(math.exp(-x.lo), np.inf)) if x.lo >= 0 else float("nan")
    e_mx = Interval(e_lo, e_hi)
    numer = one - e_mx
    tau = numer / a

    eps = (IR - IB) / IK2
    U = IUL + IW * eps.sqrt()
    E_H_cap = tau.sqr() * U.sqr() * IR
    holds = bool(E_H_cap.hi <= eps.lo + 1e-12)
    return eps, E_H_cap, tau, holds


def build_l0011_short_certificate(
    n: int = 16,
    k_ic: int = 4,
    E0: float = 0.5,
    nu: float = 0.1,
    R: float = 48.15928260103266,
    safety: float = 0.99,
) -> L0011ShortCertificate:
    K_ic2 = linf_shell_K2(k_ic)
    K2_full, M_full = full_dealias_exact_stats(n)
    M_L = shell_mode_count(n, k_ic, "linf")
    M_H = M_full - M_L
    kappa = float(kappa_H_min(n, k_ic, "linf") or 0.0)
    K_full = math.sqrt(float(K2_full))
    B = float(K_ic2) * E0
    U_L = embedding_U(M_L, E0, True)
    W = embedding_W(M_H, True)

    T_star = max_T_bootstrap(B, K_full, U_L, W, R, nu, kappa)
    T_cert = safety * T_star
    eps, E_H_cap, tau, holds = _bootstrap_interval(
        K_ic2, K2_full, M_L, M_H, E0, nu, kappa, R, T_cert
    )
    if not holds:
        raise RuntimeError("L-0011 bootstrap failed at T_cert")

    payload = json.dumps(
        {"n": n, "k_ic": k_ic, "R": R, "T_cert": T_cert, "div_free": True, "viscous": True},
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    return L0011ShortCertificate(
        certificate_id="CERT-L0011-duhamel-CS0002-short-N16",
        lemma_ref="L-0011 div-free + viscous Duhamel bootstrap",
        route="B",
        evidence_level="N5",
        n=n,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        R=R,
        T_cert=T_cert,
        T_star=T_star,
        tau_hi=float(tau.hi),
        eps_lo=float(eps.lo),
        E_H_cap_hi=float(E_H_cap.hi),
        bootstrap_holds=holds,
        clay_implication="None. Short-time finite Galerkin only; not C-S-0002 at T=0.02.",
        payload_sha256=digest,
        notes=(
            f"T_cert={T_cert:.8f} (0.99*T_star={T_star:.8f}). "
            f"tau_hi={float(tau.hi):.8e}, E_H_cap_hi={float(E_H_cap.hi):.8e} <= eps_lo={float(eps.lo):.8e}."
        ),
    )


def verify_l0011_short_certificate(cert: dict) -> tuple[bool, dict]:
    K_ic2 = linf_shell_K2(int(cert["k_ic"]))
    K2_full, M_full = full_dealias_exact_stats(int(cert["n"]))
    M_L = shell_mode_count(int(cert["n"]), int(cert["k_ic"]), "linf")
    M_H = M_full - M_L
    kappa = float(kappa_H_min(int(cert["n"]), int(cert["k_ic"]), "linf") or 0.0)
    eps, E_H_cap, tau, holds = _bootstrap_interval(
        K_ic2,
        K2_full,
        M_L,
        M_H,
        float(cert["E0"]),
        float(cert["nu"]),
        kappa,
        float(cert["R"]),
        float(cert["T_cert"]),
    )
    checks = {
        "K_ic_matches": K_ic2 == int(cert["K_ic_squared"]),
        "M_L_matches": M_L == int(cert["M_L"]),
        "M_H_matches": M_H == int(cert["M_H"]),
        "bootstrap_holds": holds,
        "stored_flag": bool(cert.get("bootstrap_holds")),
        "E_H_cap_valid_upper": float(cert["E_H_cap_hi"]) >= float(E_H_cap.hi) - 1e-9,
        "short_horizon": float(cert["T_cert"]) < 0.02,
        "evidence_level_N5": cert.get("evidence_level") == "N5",
        "no_clay_claim": "None" in cert.get("clay_implication", ""),
    }
    return all(checks.values()), checks


def save_l0011_short_certificate(cert: L0011ShortCertificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
