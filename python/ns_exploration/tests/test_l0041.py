"""Tests for L-0041 shell-6 Sym and C-R-0007."""

from __future__ import annotations

from ns_exploration.conjectures.l0041_shell6_sym import (
    SUPPORT_CR0007,
    build_cr0007,
    lemma_l0041,
)
from ns_exploration.validation.l0041_certificate import (
    build_l0041_certificate,
    verify_l0041_certificate,
)


def test_l0041_techo_and_cr0007():
    b = lemma_l0041(n=24)
    assert b.band_123456_C > b.C_dagger
    assert b.seed_12346_C <= b.C_dagger
    assert b.support_cr0007_C <= b.C_dagger
    assert 6 in (b.support_cr0007 or [])
    assert b.closes_cr0007
    assert b.consecutive_6_techo
    assert not b.closes_c0007_all_ic


def test_cr0007_build():
    cr = build_cr0007(n=24)
    assert cr["id"] == "C-R-0007"
    assert cr["shells"] == list(SUPPORT_CR0007)
    assert 6 in cr["shells"]
    assert cr["C_shor_sym"] <= cr["C_dagger"]


def test_l0041_certificate():
    cert = build_l0041_certificate(n=24)
    ok, checks = verify_l0041_certificate(cert.as_dict())
    assert ok and all(checks.values())
