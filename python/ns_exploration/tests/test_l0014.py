"""Tests for L-0014 shell→high triad short-time bound."""

from __future__ import annotations

from ns_exploration.conjectures.l0012_cancellation import omega_cancellation
from ns_exploration.conjectures.l0013_embedding import lemma_l0013
from ns_exploration.conjectures.l0014_triad import (
    effective_U_from_triads,
    lemma_l0014,
    shell_high_triad_Rmax,
)
from ns_exploration.validation.l0014_certificate import (
    build_l0014_short_certificate,
    verify_l0014_short_certificate,
)


def test_l0014_sharper_than_l0013():
    b13 = lemma_l0013()
    b14 = lemma_l0014()
    assert b14.R_H == 324
    assert b14.T_star > b13.T_star
    assert b14.U_eff < b13.U_L
    assert b14.T_star < 0.02
    assert b14.Omega_at_T_star <= b14.cs0002_M + 1e-6


def test_l0014_triad_count_and_U_eff():
    R_H = shell_high_triad_Rmax(16, 4, "linf")
    assert R_H == 324
    assert abs(effective_U_from_triads(R_H, 0.5) - 18.0) < 1e-12


def test_l0014_certificate():
    cert = build_l0014_short_certificate()
    assert cert.closes
    assert cert.R_H == 324
    assert cert.Omega_bound_hi <= cert.R + 1e-6
    ok, checks = verify_l0014_short_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0014_cannot_reach_parent_horizon():
    b = lemma_l0014()
    assert (
        omega_cancellation(0.02, b.B, b.K_full_squared**0.5, b.U_eff, b.W, b.nu, b.kappa_H)
        > b.cs0002_M
    )
    cert = build_l0014_short_certificate().as_dict()
    cert["Omega_bound_hi"] = cert["R"] + 10.0
    ok, checks = verify_l0014_short_certificate(cert)
    assert not ok
    assert not checks["Omega_le_R"]
