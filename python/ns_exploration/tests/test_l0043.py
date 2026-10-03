"""Tests for L-0043 one-pol Sym and C-R-0008."""

from __future__ import annotations

from ns_exploration.conjectures.l0043_onepol_shellblock import (
    SUPPORT_CR0008,
    cr0008_record,
    lemma_l0043,
    sym_onepol_C_shor,
)
from ns_exploration.validation.l0043_certificate import (
    build_l0043_certificate,
    verify_l0043_certificate,
)


def test_l0043_onepol_closes_cr0008():
    b = lemma_l0043(n=24)
    assert b.support_cr0008 == list(SUPPORT_CR0008)
    assert b.support_cr0008_C <= b.C_dagger
    assert b.onepol_123456_C <= b.C_dagger
    assert b.onepol_techo_12345689_C > b.C_dagger
    assert b.C_shellblock_lo > b.C_dagger
    assert b.C_kyfan_partial > b.C_dagger
    assert b.closes_cr0008
    assert not b.closes_123456_twopol
    assert not b.closes_c0007_all_ic


def test_l0043_cr0008_record():
    b = lemma_l0043(n=24)
    cr = cr0008_record(b)
    assert cr["id"] == "C-R-0008"
    assert cr["shells"] == list(SUPPORT_CR0008)
    assert cr["C_shor_sym"] <= cr["C_dagger"]


def test_l0043_sym_onepol_api():
    r = sym_onepol_C_shor(24, SUPPORT_CR0008, branch="first")
    assert r["D"] == 92
    assert r["C_shor_sym"] < 9.0


def test_l0043_certificate():
    cert = build_l0043_certificate(n=24)
    ok, checks = verify_l0043_certificate(cert.as_dict())
    assert ok and all(checks.values())
