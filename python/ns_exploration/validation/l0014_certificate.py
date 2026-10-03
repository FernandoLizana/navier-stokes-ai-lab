"""
N5 certificate for L-0014: shell→high triad multiplicity + L-0012 ODE.
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import linf_shell_K2
from ns_exploration.conjectures.l0010_duhamel_short import kappa_H_min
from ns_exploration.conjectures.l0012_cancellation import max_T_cancellation
from ns_exploration.conjectures.l0013_embedding import (
    embedding_W_vector,
    nonzero_high_mode_count,
    nonzero_shell_mode_count,
)
from ns_exploration.conjectures.l0014_triad import (
    effective_U_from_triads,
    shell_high_triad_Rmax,
    triad_C_delta,
)
from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats


@dataclass
class L0014ShortCertificate:
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
    R_H: int
    C_delta: float
    kappa_H: float
    E0: float
    nu: float
    R: float
    T_cert: float
    T_star: float
    U_eff_hi: float
    W_hi: float
    alpha_hi: float
    beta_hi: float
    z_hi: float
    Omega_bound_hi: float
    closes: bool
    clay_implication: str
    payload_sha256: str
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _exp_outward(x: Interval) -> Interval:
    return Interval(
        float(np.nextafter(math.exp(x.lo), -np.inf)),
        float(np.nextafter(math.exp(x.hi), np.inf)),
    )


def build_l0014_short_certificate(
    n: int = 16,
    k_ic: int = 4,
    E0: float = 0.5,
    nu: float = 0.1,
    R: float = 48.15928260103266,
    safety: float = 0.99,
) -> L0014ShortCertificate:
    K_ic2 = linf_shell_K2(k_ic)
    K2_full, _ = full_dealias_exact_stats(n)
    M_L = nonzero_shell_mode_count(n, k_ic, "linf")
    M_H = nonzero_high_mode_count(n, k_ic, "linf")
    R_H = shell_high_triad_Rmax(n, k_ic, "linf")
    C_delta = triad_C_delta(R_H)
    kappa = float(kappa_H_min(n, k_ic, "linf") or 0.0)
    K_full = math.sqrt(float(K2_full))

    IE0 = Interval.from_float(E0)
    Inu = Interval.from_float(nu)
    Ikap = Interval.from_float(kappa)
    IB = Interval(float(K_ic2), float(K_ic2)) * IE0
    IK2 = Interval(float(K2_full), float(K2_full))
    # U_eff = √R_H √(2E0)
    IUL = Interval(float(R_H), float(R_H)).sqrt() * (Interval(2.0, 2.0) * IE0).sqrt()
    IW = (Interval(2.0, 2.0) * Interval(float(M_H), float(M_H))).sqrt()

    sqrtB = IB.sqrt()
    beta = IUL * sqrtB
    alpha = IW * sqrtB - Inu * Ikap

    B = float(K_ic2) * E0
    U_eff = effective_U_from_triads(R_H, E0)
    Wf = embedding_W_vector(M_H)
    T_star = max_T_cancellation(B, K_full, U_eff, Wf, R, nu, kappa)
    T_cert = safety * T_star
    IT = Interval.from_float(T_cert)

    if alpha.lo <= 0:
        raise RuntimeError("expected alpha > 0 for C-S N=16 regime")
    exp_at = _exp_outward(alpha * IT)
    z = (beta / alpha) * (exp_at - Interval(1.0, 1.0))
    Omega = IB + IK2 * z.sqr()

    closes = bool(Omega.hi <= R + 1e-9)
    if not closes:
        raise RuntimeError(f"Omega.hi={Omega.hi} exceeds R={R} at T_cert")

    payload = json.dumps(
        {
            "n": n,
            "k_ic": k_ic,
            "R": R,
            "T_cert": T_cert,
            "R_H": R_H,
            "lemma": "L-0014",
        },
        sort_keys=True,
    ).encode()
    digest = hashlib.sha256(payload).hexdigest()

    return L0014ShortCertificate(
        certificate_id="CERT-L0014-triad-CS0002-short-N16",
        lemma_ref="L-0014 shell→high triad multiplicity + L-0012 ODE",
        route="B",
        evidence_level="N5",
        n=n,
        k_ic=k_ic,
        K_ic_squared=K_ic2,
        K_full_squared=K2_full,
        M_L=M_L,
        M_H=M_H,
        R_H=R_H,
        C_delta=C_delta,
        kappa_H=kappa,
        E0=E0,
        nu=nu,
        R=R,
        T_cert=T_cert,
        T_star=T_star,
        U_eff_hi=float(IUL.hi),
        W_hi=float(IW.hi),
        alpha_hi=float(alpha.hi),
        beta_hi=float(beta.hi),
        z_hi=float(z.hi),
        Omega_bound_hi=float(Omega.hi),
        closes=closes,
        clay_implication="None. Short-time finite Galerkin only; not C-S-0002 at T=0.02.",
        payload_sha256=digest,
        notes=(
            f"R_H={R_H}; C_△={C_delta:.8f}; T_cert={T_cert:.8f}; "
            f"Omega_hi={float(Omega.hi):.8f} <= R={R:.8f}."
        ),
    )


def verify_l0014_short_certificate(cert: dict) -> tuple[bool, dict]:
    rebuilt = build_l0014_short_certificate(
        n=int(cert["n"]),
        k_ic=int(cert["k_ic"]),
        E0=float(cert["E0"]),
        nu=float(cert["nu"]),
        R=float(cert["R"]),
        safety=float(cert["T_cert"]) / float(cert["T_star"])
        if float(cert["T_star"]) > 0
        else 0.99,
    )
    checks = {
        "K_ic_matches": rebuilt.K_ic_squared == int(cert["K_ic_squared"]),
        "M_L_matches": rebuilt.M_L == int(cert["M_L"]),
        "M_H_matches": rebuilt.M_H == int(cert["M_H"]),
        "R_H_matches": rebuilt.R_H == int(cert["R_H"]),
        "closes": rebuilt.closes and bool(cert.get("closes")),
        "Omega_valid_upper": float(cert["Omega_bound_hi"]) >= rebuilt.Omega_bound_hi - 1e-6,
        "Omega_le_R": float(cert["Omega_bound_hi"]) <= float(cert["R"]) + 1e-6,
        "short_horizon": float(cert["T_cert"]) < 0.02,
        "evidence_level_N5": cert.get("evidence_level") == "N5",
        "no_clay_claim": "None" in cert.get("clay_implication", ""),
        "sharper_than_L0013_T": float(cert["T_star"]) > 0.0032350612481493943,
    }
    return all(checks.values()), checks


def save_l0014_short_certificate(cert: L0014ShortCertificate, path: str | Path) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(cert.as_dict(), indent=2), encoding="utf-8")
