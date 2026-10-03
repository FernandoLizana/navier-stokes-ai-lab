"""Tests for L-0022 gamma bootstrap."""

from __future__ import annotations

from ns_exploration.conjectures.l0020_duhamel_h1 import lemma_l0020
from ns_exploration.conjectures.l0022_gamma_bootstrap import gamma_crit, lemma_l0022
from ns_exploration.validation.l0022_certificate import (
    build_l0022_certificate,
    verify_l0022_certificate,
)


def test_l0022_gamma_crit_arithmetic():
    b20 = lemma_l0020()
    g = gamma_crit(0.5, b20.c0005_M, b20.Ncrit_c0005)
    assert abs(g - 1.746882907641106) < 1e-8
    b = lemma_l0022(empirical=False)
    assert abs(b.gamma_crit - g) < 1e-12
    assert not b.proved_closes_c0005


def test_l0022_certificate():
    cert = build_l0022_certificate(empirical=False)
    ok, checks = verify_l0022_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0022_adversarial_N_exceeds_Ncrit():
    """Uniform N_* ≤ N_crit is false for all-IC (VJP ascent)."""
    b = lemma_l0022(empirical=True, n_random=3, ascent_steps=12)
    assert b.N_emp_max > b.Ncrit_c0005
    assert b.gamma_emp_max > b.gamma_crit
    assert not b.uniform_hypothesis_holds_empirically
    # Still far from refuting C-0005 on these fields
    assert b.Omega_T_emp_max < b.c0005_M
