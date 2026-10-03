"""
Tests for the L-0003 uniform interval certificate (shell + full dealias).
"""

from __future__ import annotations

from ns_exploration.validation.l0003_certificate import (
    build_l0003_certificate,
    build_l0003_full_dealias_certificate,
    full_dealias_exact_stats,
)
from ns_exploration.validation.l0003_certificate_verify import verify_l0003_certificate
from ns_exploration.validation.l0006_certificate import build_certificate, shell_exact_stats


def test_shell_k3_stats():
    K2, M = shell_exact_stats(12, 3)
    assert K2 == 9
    assert M == 123


def test_full_dealias_n12_stats():
    K2, M = full_dealias_exact_stats(12)
    assert K2 == 27  # (3,3,3)
    assert M == 343  # 7^3 with |k_i| < 4


def test_l0006_algebraic_fails_at_k3():
    """Documents why we fall back to L-0003 for k_shell=3."""
    cert = build_certificate(n=12, k_shell=3, E0=0.5, t=0.02)
    assert not cert.closes
    assert cert.Omega_bound_hi == float("inf")


def test_l0003_certificate_closes_and_bounds():
    cert = build_l0003_certificate(n=12, k_shell=3, E0=0.5, nu=0.1)
    assert cert.closes
    assert cert.support == "shell"
    assert cert.K_squared_int == 9
    assert cert.M == 123
    # Ω0 = 4.5, Ω_eq = 18450 — bound is the equilibrium
    assert abs(cert.Omega0_hi - 4.5) < 1e-9
    assert 18449.0 < cert.Omega_bound_hi < 18451.0
    assert cert.evidence_level == "N5"
    assert "None" in cert.clay_implication


def test_l0003_full_dealias_certificate():
    cert = build_l0003_full_dealias_certificate(n=12, E0=0.5, nu=0.1)
    assert cert.closes
    assert cert.support == "full_dealias"
    assert cert.k_shell is None
    assert cert.K_squared_int == 27
    assert cert.M == 343
    assert abs(cert.Omega0_hi - 13.5) < 1e-9
    assert 51449.0 < cert.Omega_bound_hi < 51451.0
    ok, checks = verify_l0003_certificate(cert.as_dict())
    assert ok
    assert all(checks.values())


def test_l0003_verifier_accepts():
    ok, checks = verify_l0003_certificate(build_l0003_certificate().as_dict())
    assert ok
    assert all(checks.values())


def test_l0003_verifier_rejects_understated_bound():
    cert = build_l0003_certificate().as_dict()
    cert["Omega_bound_hi"] = 1.0
    ok, checks = verify_l0003_certificate(cert)
    assert not ok
    assert not checks["bound_is_valid_upper"]


def test_l0003_verifier_rejects_wrong_M():
    cert = build_l0003_certificate().as_dict()
    cert["M"] = 1
    ok, checks = verify_l0003_certificate(cert)
    assert not ok
    assert not checks["M_matches"]


def test_l0003_full_verifier_rejects_shell_stats():
    """Full-dealias cert must not accept shell (K2,M)."""
    cert = build_l0003_full_dealias_certificate().as_dict()
    cert["K_squared_int"] = 9
    cert["M"] = 123
    ok, checks = verify_l0003_certificate(cert)
    assert not ok
    assert not checks["K_squared_matches"] or not checks["M_matches"]
