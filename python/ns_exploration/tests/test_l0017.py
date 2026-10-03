"""Tests for L-0017 radial H×L cross short-time bound."""

from __future__ import annotations

from ns_exploration.conjectures.l0016_radial import lemma_l0016
from ns_exploration.conjectures.l0017_cross import (
    gamma_cross_from_rho,
    lemma_l0017,
    omega_cancellation_cross,
)
from ns_exploration.validation.l0017_certificate import (
    build_l0017_short_certificate,
    verify_l0017_short_certificate,
)


def test_l0017_sharper_than_l0016():
    b16 = lemma_l0016()
    b17 = lemma_l0017()
    assert b17.gamma_cross < b17.W_sqrt_B
    assert b17.T_star > b16.T_star
    assert b17.T_star < 0.02
    assert b17.Omega_at_T_star <= b17.cs0002_M + 1e-6
    assert abs(b17.gamma_cross - gamma_cross_from_rho(1392.0, 0.5)) < 1e-12


def test_l0017_certificate():
    cert = build_l0017_short_certificate()
    assert cert.closes
    assert abs(cert.rho_star - 1392.0) < 1e-9
    assert cert.Omega_bound_hi <= cert.R + 1e-6
    ok, checks = verify_l0017_short_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0017_cannot_reach_parent_horizon():
    b = lemma_l0017()
    assert (
        omega_cancellation_cross(
            0.02, b.B, b.K_full_squared**0.5, b.U_eff, b.gamma_cross, b.W, b.nu, b.kappa_H
        )
        > b.cs0002_M
    )
    cert = build_l0017_short_certificate().as_dict()
    cert["Omega_bound_hi"] = cert["R"] + 10.0
    ok, checks = verify_l0017_short_certificate(cert)
    assert not ok
    assert not checks["Omega_le_R"]


def test_l0017_near_parent_horizon():
    b = lemma_l0017()
    assert b.T_star > 0.015
    assert b.T_star / 0.02 > 0.75
