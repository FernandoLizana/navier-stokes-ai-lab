"""
Independent minimal verifier for CERT-L0003-* (shell or full_dealias).

Re-derives (K², M) from first principles according to `support`, recomputes the
outward enclosure of max(K² E0, 6 M E0² / ν²), and checks that the stored
Omega_bound_hi is a valid upper bound.
"""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0003_certificate import full_dealias_exact_stats
from ns_exploration.validation.l0006_certificate import shell_exact_stats


def verify_l0003_certificate(cert: dict) -> tuple[bool, dict]:
    checks: dict[str, bool] = {}
    support = cert.get("support", "shell")  # backward compatible with older JSONs

    if support == "shell":
        k_shell = cert.get("k_shell")
        if k_shell is None:
            return False, {"support_ok": False}
        K2, M = shell_exact_stats(int(cert["n"]), int(k_shell))
    elif support == "full_dealias":
        K2, M = full_dealias_exact_stats(int(cert["n"]))
    else:
        return False, {"support_ok": False}

    checks["support_ok"] = True
    checks["K_squared_matches"] = K2 == int(cert["K_squared_int"])
    checks["M_matches"] = M == int(cert["M"])

    IE0 = Interval.from_float(float(cert["E0"]))
    Inu = Interval.from_float(float(cert["nu"]))
    Omega0 = Interval(float(K2), float(K2)) * IE0
    Omega_eq = (Interval(6.0, 6.0) * Interval(float(M), float(M)) * IE0.sqr()) / Inu.sqr()

    recomputed_hi = max(float(Omega0.hi), float(Omega_eq.hi))
    stored = float(cert["Omega_bound_hi"])

    checks["Omega0_hi_consistent"] = float(cert["Omega0_hi"]) >= float(Omega0.hi) - 1e-12
    checks["Omega_eq_hi_consistent"] = float(cert["Omega_eq_hi"]) >= float(Omega_eq.hi) - 1e-12
    checks["bound_is_valid_upper"] = stored >= recomputed_hi - 1e-12
    checks["closes"] = bool(cert.get("closes", False))
    checks["evidence_level_N5"] = cert.get("evidence_level") == "N5"
    checks["no_clay_claim"] = "None" in cert.get("clay_implication", "")

    return all(checks.values()), checks


def verify_l0003_file(path: str | Path) -> tuple[bool, dict]:
    cert = json.loads(Path(path).read_text(encoding="utf-8"))
    return verify_l0003_certificate(cert)
