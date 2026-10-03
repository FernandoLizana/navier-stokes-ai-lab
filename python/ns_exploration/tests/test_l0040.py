"""Tests for L-0040 Sym Shor and C-R-0006."""

from __future__ import annotations

from ns_exploration.conjectures.l0040_sym_shor import (
    SUPPORT_CR0006,
    build_cr0006,
    lemma_l0040,
    sym_band_C_shor,
)
from ns_exploration.validation.l0040_certificate import (
    build_l0040_certificate,
    verify_l0040_certificate,
)


def test_l0040_sym_closes_12345_not_123456():
    b = lemma_l0040(n=24)
    assert b.band_12345_C <= b.C_dagger
    assert b.band_123456_C > b.C_dagger
    assert b.band_1234_C_sym < b.band_12345_C
    assert b.closes_cr0006
    assert not b.closes_c0007_all_ic


def test_l0040_sym_values():
    d = sym_band_C_shor(24, SUPPORT_CR0006)
    assert d["D"] == 112
    assert 9.0 < d["C_shor_sym"] < 9.1


def test_cr0006_build():
    cr = build_cr0006(n=24)
    assert cr["id"] == "C-R-0006"
    assert cr["shells"] == [1, 2, 3, 4, 5]
    assert cr["C_shor_sym"] <= cr["C_dagger"]


def test_l0040_certificate():
    cert = build_l0040_certificate(n=24)
    ok, checks = verify_l0040_certificate(cert.as_dict())
    assert ok and all(checks.values())
