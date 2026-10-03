"""Tests for L-0011 sharpened short-time bounds."""

from __future__ import annotations

from ns_exploration.conjectures.l0010_duhamel_short import lemma_l0010
from ns_exploration.conjectures.l0011_sharpened_short import lemma_l0011
from ns_exploration.validation.l0011_certificate import (
    build_l0011_short_certificate,
    verify_l0011_short_certificate,
)


def test_l0011_improves_Tstar():
    b0 = lemma_l0010()
    b1 = lemma_l0011()
    assert b1.T_star_rigorous > b0.T_duhamel_for_cs0002
    assert b1.T_star_rigorous < 0.02
    assert 1.15 < b1.improvement_vs_L0010 < 1.5


def test_l0011_LxL_conditional_longer():
    b1 = lemma_l0011()
    assert b1.T_star_LxL_conditional > b1.T_star_rigorous
    assert b1.Omega_LxL_at_T_target > b1.cs0002_M  # still fails at 0.02
    assert b1.LxL_evidence == "N6"


def test_l0011_certificate():
    cert = build_l0011_short_certificate()
    assert cert.bootstrap_holds
    ok, checks = verify_l0011_short_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0011_rejects_long_T():
    cert = build_l0011_short_certificate().as_dict()
    cert["T_cert"] = 0.02
    ok, checks = verify_l0011_short_certificate(cert)
    assert not ok
    assert not checks["bootstrap_holds"]
