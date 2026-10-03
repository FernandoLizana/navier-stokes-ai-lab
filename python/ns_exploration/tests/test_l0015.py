"""Tests for L-0015 S_max short-time bound."""

from __future__ import annotations

from ns_exploration.conjectures.l0012_cancellation import omega_cancellation
from ns_exploration.conjectures.l0014_triad import lemma_l0014
from ns_exploration.conjectures.l0015_smax import lemma_l0015, shell_high_Smax
from ns_exploration.validation.l0015_certificate import (
    build_l0015_short_certificate,
    verify_l0015_short_certificate,
)


def test_l0015_sharper_than_l0014():
    b14 = lemma_l0014()
    b15 = lemma_l0015()
    assert abs(b15.S_max - 6132.0) < 1e-9
    assert b15.T_star > b14.T_star
    assert b15.U_eff < b14.U_eff
    assert b15.T_star < 0.02
    assert b15.Omega_at_T_star <= b15.cs0002_M + 1e-6


def test_l0015_smax_value():
    assert abs(shell_high_Smax(16, 4, "linf") - 6132.0) < 1e-9


def test_l0015_certificate():
    cert = build_l0015_short_certificate()
    assert cert.closes
    assert abs(cert.S_max - 6132.0) < 1e-9
    assert cert.Omega_bound_hi <= cert.R + 1e-6
    ok, checks = verify_l0015_short_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0015_cannot_reach_parent_horizon():
    b = lemma_l0015()
    assert (
        omega_cancellation(0.02, b.B, b.K_full_squared**0.5, b.U_eff, b.W, b.nu, b.kappa_H)
        > b.cs0002_M
    )
    cert = build_l0015_short_certificate().as_dict()
    cert["Omega_bound_hi"] = cert["R"] + 10.0
    ok, checks = verify_l0015_short_certificate(cert)
    assert not ok
    assert not checks["Omega_le_R"]
