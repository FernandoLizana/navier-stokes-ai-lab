"""Tests for L-0038 band SVD Shor and C-R-0004."""

from __future__ import annotations

from ns_exploration.conjectures.l0038_band_svd_shor import (
    build_cr0004,
    exact_band_C_shor,
    lemma_l0038,
)
from ns_exploration.validation.l0038_certificate import (
    build_l0038_certificate,
    verify_l0038_certificate,
)


def test_l0038_bands_meet_and_fail():
    b = lemma_l0038(n=24)
    assert b.band_123_C <= b.C_dagger
    assert b.band_12_C <= b.C_dagger
    assert b.band_125_C <= b.C_dagger
    assert b.band_1234_C > b.C_dagger
    assert b.closes_cr0004
    assert not b.closes_c0007_all_ic


def test_l0038_exact_svd_stable():
    d = exact_band_C_shor(24, (1, 2, 3))
    assert d["D"] == 52
    assert d["C_shor"] < 9.2


def test_cr0004_build():
    cr = build_cr0004(n=24)
    assert cr["id"] == "C-R-0004"
    assert cr["shells"] == [1, 2, 3]
    assert cr["C_shor"] <= cr["C_dagger"]


def test_l0038_certificate():
    cert = build_l0038_certificate(n=24)
    ok, checks = verify_l0038_certificate(cert.as_dict())
    assert ok and all(checks.values())
