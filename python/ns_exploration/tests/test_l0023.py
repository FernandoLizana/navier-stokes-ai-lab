"""Tests for L-0023 cubic-dissipation bootstrap."""

from __future__ import annotations

from ns_exploration.conjectures.l0023_cubic_dissipation import (
    C_star_for_M,
    lemma_l0023,
    omega_ode_bound,
)
from ns_exploration.validation.l0023_certificate import (
    build_l0023_certificate,
    verify_l0023_certificate,
)


def test_l0023_cstar_closes_exactly_at_threshold():
    E0, nu, T = 0.5, 0.1, 0.02
    Omega0 = 147 * E0
    M = 61.23693461895651
    Cstar = C_star_for_M(E0, nu, T, Omega0, M)
    assert 2.0 < Cstar < 2.3
    om = omega_ode_bound(Cstar, E0, nu, T, Omega0)
    assert abs(om - M) < 1e-4
    # slightly larger C exceeds M
    assert omega_ode_bound(Cstar + 0.05, E0, nu, T, Omega0) > M


def test_l0023_certificate():
    cert = build_l0023_certificate(empirical=False)
    ok, checks = verify_l0023_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0023_empirics_beat_cstar():
    b = lemma_l0023(empirical=True, n_random=6, ascent_steps=8)
    assert b.C_emp_max < 0.5
    assert b.C_emp_max < b.C_star
    assert b.emp_closes
    assert not b.proved_closes_c0005
