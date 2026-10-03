"""Tests for L-0012 NS-cancellation short-time bound."""

from __future__ import annotations

from ns_exploration.conjectures.l0011_sharpened_short import lemma_l0011
from ns_exploration.conjectures.l0012_cancellation import lemma_l0012, omega_cancellation
from ns_exploration.validation.l0012_certificate import (
    build_l0012_short_certificate,
    verify_l0012_short_certificate,
)


def test_l0012_improves_Tstar():
    b11 = lemma_l0011()
    b12 = lemma_l0012()
    assert b12.T_star > b11.T_star_rigorous
    assert b12.T_star < 0.02
    assert b12.improvement_vs_L0011 > 1.5
    assert b12.Omega_at_T_star <= b12.cs0002_M + 1e-6


def test_l0012_formula_at_zero():
    b = lemma_l0012()
    assert abs(omega_cancellation(0.0, b.B, b.K_full_squared**0.5, b.U_L, b.W, b.nu, b.kappa_H) - b.B) < 1e-12


def test_l0012_certificate():
    cert = build_l0012_short_certificate()
    assert cert.closes
    assert cert.Omega_bound_hi <= cert.R + 1e-6
    ok, checks = verify_l0012_short_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0012_cannot_reach_parent_horizon():
    """Parent C-S-0002 asks T=0.02; L-0012 majorant exceeds R there."""
    b = lemma_l0012()
    assert omega_cancellation(0.02, b.B, b.K_full_squared**0.5, b.U_L, b.W, b.nu, b.kappa_H) > b.cs0002_M
    cert = build_l0012_short_certificate().as_dict()
    cert["Omega_bound_hi"] = cert["R"] + 10.0
    ok, checks = verify_l0012_short_certificate(cert)
    assert not ok
    assert not checks["Omega_le_R"]