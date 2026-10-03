"""Tests for L-0020 Stokes–Duhamel H¹ majorant."""

from __future__ import annotations

from ns_exploration.conjectures.l0019_stokes_majorant import lemma_l0019
from ns_exploration.conjectures.l0020_duhamel_h1 import (
    Ncrit_for_M,
    lemma_l0020,
    omega_duhamel_H1,
)
from ns_exploration.validation.l0020_certificate import (
    build_l0020_certificate,
    verify_l0020_certificate,
)


def test_l0020_matches_stokes_floor():
    b19 = lemma_l0019()
    b20 = lemma_l0020()
    assert abs(b20.Stokes_floor - b19.Stokes_floor) < 1e-10
    assert b20.rho_star == 6192.0
    assert b20.r_star == 86


def test_l0020_young_does_not_close_c0005():
    b = lemma_l0020()
    assert not b.young_closes_c0005
    assert b.Omega_young_H1 > b.envelope_cap  # worse than L-0018
    assert b.Ncrit_c0005 > 9.0
    assert b.Ncrit_c0005 < 10.0


def test_l0020_ncrit_formula():
    b = lemma_l0020()
    # At N=Ncrit, bound hits M
    om = omega_duhamel_H1(b.S_T, b.E0, b.Ncrit_c0005, b.I_sigma)
    assert abs(om - b.c0005_M) < 1e-6
    assert abs(b.Ncrit_c0005 - Ncrit_for_M(b.S_T, b.E0, b.I_sigma, b.c0005_M)) < 1e-12


def test_l0020_certificate():
    cert = build_l0020_certificate()
    ok, checks = verify_l0020_certificate(cert.as_dict())
    assert ok and all(checks.values())
    # tamper
    d = cert.as_dict()
    d["Ncrit_c0005"] = 0.0
    ok2, ch2 = verify_l0020_certificate(d)
    assert not ok2 and not ch2["Ncrit"]
