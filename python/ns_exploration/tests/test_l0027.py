"""Tests for L-0027 / C-R-0002."""

from __future__ import annotations

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.validation.l0027_certificate import (
    build_l0027_certificate,
    verify_l0027_certificate,
)


def test_l0027_low_slab_conditional():
    b = lemma_l0027(n=24, empirical=False, n_grid=21)
    assert b.C_dagger > 5.0
    assert b.Omega_T_at_Cdagger <= b.c0007_M + 1e-9
    assert b.all_ic_cubic_C0_floor > b.c0007_M
    assert not b.all_ic_cubic_closes_c0007
    assert not b.proved_closes_c0007


def test_l0027_certificate():
    cert = build_l0027_certificate(n=24)
    ok, checks = verify_l0027_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_cr0002_high_slab_matches_l0026():
    b = lemma_l0026(n=24, n_grid=41)
    assert b.high_slab_closes
    assert b.Omega_T_at_star <= b.c0007_M + 1e-9
    assert not b.closes_c0007
