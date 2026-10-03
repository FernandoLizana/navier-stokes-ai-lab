"""Tests for L-0046 greedy physical fullsym and C-R-0011 (slow: ~15–20 min)."""

from __future__ import annotations

from ns_exploration.conjectures.l0046_greedy_fullsym import (
    SUPPORT_CR0011,
    cr0011_record,
    lemma_l0046,
)
from ns_exploration.validation.l0046_certificate import (
    build_l0046_certificate,
    verify_l0046_certificate,
)


def test_l0046_cr0011_and_certificate():
    b = lemma_l0046(n=24)
    assert b.support_cr0011 == list(SUPPORT_CR0011)
    assert 147 in b.support_cr0011
    assert b.C_fullsym_cr0011 <= b.C_dagger
    assert b.closes_cr0011
    assert not b.closes_c0007_all_ic
    cr = cr0011_record(b)
    assert cr["id"] == "C-R-0011"
    cert = build_l0046_certificate(n=24, bound=b)
    ok, checks = verify_l0046_certificate(cert.as_dict())
    assert ok and all(checks.values())
