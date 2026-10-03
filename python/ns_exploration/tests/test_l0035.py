"""Tests for L-0035 shell two-point ‖∇u‖_∞ stretch techo."""

from __future__ import annotations

from ns_exploration.conjectures.l0035_ginf_stretch import lemma_l0035
from ns_exploration.validation.l0035_certificate import (
    build_l0035_certificate,
    verify_l0035_certificate,
)


def test_l0035_meets_cdagger_only_at_equipartition():
    b = lemma_l0035(n=24, n_grid=7)
    assert b.C_at_equipartition <= b.C_dagger
    assert b.C_max_low_slab > b.C_dagger
    assert abs(b.Omega_c - 0.5) < 0.05
    assert b.ginf_ode_worst > b.c0007_M
    assert b.minprod_ode_worst > b.c0007_M
    assert not b.closes_c0007


def test_l0035_certificate():
    cert = build_l0035_certificate(n=24)
    ok, checks = verify_l0035_certificate(cert.as_dict())
    assert ok and all(checks.values())
