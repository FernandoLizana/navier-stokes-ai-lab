"""
Independent minimal verifier for CERT-L0006-shell-k2.

Design goal (§25): the verifier is *small* and re-derives the guaranteed
upper bound from the exact integer data stored in the certificate, then
checks that the stored Omega_bound_hi is >= the freshly recomputed enclosure.
It does NOT trust the generator's numeric fields other than the exact
integers (K², M) and the rational inputs (E0, t, k_shell, n).
"""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.validation.intervals import Interval
from ns_exploration.validation.l0006_certificate import shell_exact_stats


def verify_certificate(cert: dict) -> tuple[bool, dict]:
    checks: dict[str, bool] = {}

    # 1. Re-derive exact integer data from first principles.
    K2, M = shell_exact_stats(int(cert["n"]), int(cert["k_shell"]))
    checks["K_squared_matches"] = K2 == int(cert["K_squared_int"])
    checks["M_matches"] = M == int(cert["M"])

    # 2. Independently recompute the outward enclosure.
    IE0 = Interval.from_float(float(cert["E0"]))
    It = Interval.from_float(float(cert["t"]))
    Omega0 = Interval(float(K2), float(K2)) * IE0
    a = (Interval(2.0, 2.0) * Interval(float(3 * M), float(3 * M))).sqrt()
    prod = a * Omega0.sqrt() * It

    checks["prod_hi_consistent"] = prod.hi <= float(cert["prod_hi"]) + 1e-15
    closes = prod.hi < 1.0
    checks["closes_matches"] = closes == bool(cert["closes"])

    if closes:
        denom = Interval(1.0, 1.0) - prod
        cap = Omega0 / denom.sqr()
        # Stored bound must not underclaim the true guaranteed ceiling.
        checks["bound_is_valid_upper"] = float(cert["Omega_bound_hi"]) >= cap.hi - 1e-12
    else:
        checks["bound_is_valid_upper"] = cert["Omega_bound_hi"] == float("inf")

    checks["evidence_level_N5"] = cert.get("evidence_level") == "N5"
    checks["no_clay_claim"] = "None" in cert.get("clay_implication", "")

    ok = all(checks.values())
    return ok, checks


def verify_file(path: str | Path) -> tuple[bool, dict]:
    cert = json.loads(Path(path).read_text(encoding="utf-8"))
    return verify_certificate(cert)
