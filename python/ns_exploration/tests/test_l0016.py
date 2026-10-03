"""Tests for L-0016 radial Young short-time bound."""

from __future__ import annotations

from ns_exploration.conjectures.l0012_cancellation import omega_cancellation
from ns_exploration.conjectures.l0015_smax import lemma_l0015
from ns_exploration.conjectures.l0016_radial import lemma_l0016, shell_radial_rho_star
from ns_exploration.validation.l0016_certificate import (
    build_l0016_short_certificate,
    verify_l0016_short_certificate,
)


def test_l0016_sharper_than_l0015():
    b15 = lemma_l0015()
    b16 = lemma_l0016()
    assert abs(b16.rho_star - 1392.0) < 1e-9
    assert b16.T_star > b15.T_star
    assert b16.U_eff < b15.U_eff
    assert b16.T_star < 0.02
    assert b16.Omega_at_T_star <= b16.cs0002_M + 1e-6


def test_l0016_rho_star():
    rho, r, m = shell_radial_rho_star(16, 4, "linf")
    assert abs(rho - 1392.0) < 1e-9
    assert r == 29 and m == 48


def test_l0016_certificate():
    cert = build_l0016_short_certificate()
    assert cert.closes
    assert abs(cert.rho_star - 1392.0) < 1e-9
    assert cert.Omega_bound_hi <= cert.R + 1e-6
    ok, checks = verify_l0016_short_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0016_cannot_reach_parent_horizon():
    b = lemma_l0016()
    assert (
        omega_cancellation(0.02, b.B, b.K_full_squared**0.5, b.U_eff, b.W, b.nu, b.kappa_H)
        > b.cs0002_M
    )
    cert = build_l0016_short_certificate().as_dict()
    cert["Omega_bound_hi"] = cert["R"] + 10.0
    ok, checks = verify_l0016_short_certificate(cert)
    assert not ok
    assert not checks["Omega_le_R"]
