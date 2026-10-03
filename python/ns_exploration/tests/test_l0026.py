"""Tests for L-0026 / C-0007 techo."""

from __future__ import annotations

from ns_exploration.conjectures.l0026_c0007_techo import lemma_l0026
from ns_exploration.validation.l0026_certificate import (
    build_l0026_certificate,
    verify_l0026_certificate,
)


def test_l0026_high_slab_closes_low_open():
    b = lemma_l0026(n=24, n_grid=41)
    assert b.monotone_in_Omega0
    assert b.high_slab_closes
    assert not b.closes_c0007
    assert b.Stokes_floor < b.c0007_M < b.L0024_majorant
    assert b.E0 <= b.Omega_star < b.L0024_majorant
    assert b.Omega_T_at_star <= b.c0007_M + 1e-9
    assert b.low_slab_L0024_worst > b.c0007_M


def test_l0026_certificate():
    cert = build_l0026_certificate(n=24)
    ok, checks = verify_l0026_certificate(cert.as_dict())
    assert ok and all(checks.values())
