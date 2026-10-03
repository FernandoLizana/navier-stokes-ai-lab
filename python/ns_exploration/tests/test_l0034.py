"""Tests for L-0034 shell-factorized ‖N‖ bound."""

from __future__ import annotations

from ns_exploration.conjectures.l0034_shell_N_bound import lemma_l0034
from ns_exploration.validation.l0034_certificate import (
    build_l0034_certificate,
    verify_l0034_certificate,
)


def test_l0034_beats_young_but_not_cdagger():
    b = lemma_l0034(n=24, n_grid=9)
    assert b.N_at_equipartition < b.N_young
    assert b.C_eff_at_equipartition > b.C_dagger
    assert b.low_slab_ode_worst > b.c0007_M
    assert not b.closes_c0007


def test_l0034_certificate():
    cert = build_l0034_certificate(n=24)
    ok, checks = verify_l0034_certificate(cert.as_dict())
    assert ok and all(checks.values())
