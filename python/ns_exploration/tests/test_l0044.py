"""Tests for L-0044 physical Hermitian fullsym and C-R-0009."""

from __future__ import annotations

from ns_exploration.conjectures.l0044_physical_fullsym import (
    SUPPORT_CR0009,
    cr0009_record,
    lemma_l0044,
)
from ns_exploration.validation.l0044_certificate import (
    build_l0044_certificate,
    verify_l0044_certificate,
)


def test_l0044_fullsym_closes_cr0009():
    b = lemma_l0044(n=24)
    assert b.support_cr0009 == list(SUPPORT_CR0009)
    assert b.C_fullsym_123456 <= b.C_dagger
    assert b.C_sym_123456 > b.C_dagger
    assert b.fft_tensor_err < 1e-11
    assert b.closes_cr0009
    assert not b.closes_c0007_all_ic


def test_l0044_cr0009_record():
    b = lemma_l0044(n=24)
    cr = cr0009_record(b)
    assert cr["id"] == "C-R-0009"
    assert cr["shells"] == list(SUPPORT_CR0009)
    assert cr["C_fullsym"] <= cr["C_dagger"]


def test_l0044_certificate():
    cert = build_l0044_certificate(n=24)
    ok, checks = verify_l0044_certificate(cert.as_dict())
    assert ok and all(checks.values())
