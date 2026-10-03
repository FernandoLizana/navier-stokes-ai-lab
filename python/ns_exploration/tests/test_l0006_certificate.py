"""
Tests for the L-0006 interval certificate and its independent verifier.

These are fast, deterministic, and finite-dimensional. They check:
  - exact integer shell data (K², M),
  - the outward enclosure closes and yields a finite upper bound,
  - the stored bound is a valid upper bound (verifier passes),
  - tampering with the stored bound downward is caught.
"""

from __future__ import annotations

from ns_exploration.validation.l0006_certificate import (
    build_certificate,
    shell_exact_stats,
)
from ns_exploration.validation.l0006_certificate_verify import verify_certificate


def test_shell_exact_stats_k2():
    K2, M = shell_exact_stats(12, 2)
    assert K2 == 4  # max |k|^2 with |k|<=2 euclidean is (2,0,0)->4
    assert M == 33  # integer lattice points in Euclidean ball radius 2 (incl origin)


def test_certificate_closes_and_bounds():
    cert = build_certificate(n=12, k_shell=2, E0=0.5, t=0.02)
    assert cert.closes
    assert cert.prod_hi < 1.0
    assert 4.0 < cert.Omega_bound_hi < 8.0  # ~5.52
    assert cert.evidence_level == "N5"
    assert "None" in cert.clay_implication


def test_verifier_accepts_valid_certificate():
    cert = build_certificate()
    ok, checks = verify_certificate(cert.as_dict())
    assert ok
    assert all(checks.values())


def test_verifier_rejects_understated_bound():
    cert = build_certificate().as_dict()
    cert["Omega_bound_hi"] = 1.0  # too small -> not a valid upper bound
    ok, checks = verify_certificate(cert)
    assert not ok
    assert not checks["bound_is_valid_upper"]


def test_verifier_rejects_wrong_mode_count():
    cert = build_certificate().as_dict()
    cert["M"] = 999
    ok, checks = verify_certificate(cert)
    assert not ok
    assert not checks["M_matches"]
