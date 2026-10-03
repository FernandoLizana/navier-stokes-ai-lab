"""Tests for L-0013 vector-CS embedding short-time bound."""

from __future__ import annotations

from ns_exploration.conjectures.l0012_cancellation import lemma_l0012
from ns_exploration.conjectures.l0013_embedding import (
    embedding_U_vector,
    embedding_W_vector,
    lemma_l0013,
    nonzero_shell_mode_count,
)
from ns_exploration.validation.l0013_certificate import (
    build_l0013_short_certificate,
    verify_l0013_short_certificate,
)


def test_l0013_sharper_than_l0012():
    b12 = lemma_l0012()
    b13 = lemma_l0013()
    assert b13.T_star > b12.T_star
    assert b13.U_L < b12.U_L
    assert b13.W < b12.W
    assert b13.T_star < 0.02
    assert b13.Omega_at_T_star <= b13.cs0002_M + 1e-6


def test_l0013_vector_cs_formula():
    M = nonzero_shell_mode_count(16, 4, "linf")
    assert M == 728
    assert abs(embedding_U_vector(M, 0.5) - (M**0.5)) < 1e-12
    assert abs(embedding_W_vector(602) - (2 * 602) ** 0.5) < 1e-12


def test_l0013_certificate():
    cert = build_l0013_short_certificate()
    assert cert.closes
    assert cert.Omega_bound_hi <= cert.R + 1e-6
    assert cert.M_L == 728
    ok, checks = verify_l0013_short_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0013_cannot_reach_parent_horizon():
    b = lemma_l0013()
    from ns_exploration.conjectures.l0012_cancellation import omega_cancellation

    assert (
        omega_cancellation(0.02, b.B, b.K_full_squared**0.5, b.U_L, b.W, b.nu, b.kappa_H)
        > b.cs0002_M
    )
    cert = build_l0013_short_certificate().as_dict()
    cert["Omega_bound_hi"] = cert["R"] + 10.0
    ok, checks = verify_l0013_short_certificate(cert)
    assert not ok
    assert not checks["Omega_le_R"]
