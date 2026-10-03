"""Tests for L-0047 extended greedy fullsym and C-R-0012 (slow: ~10–20 min)."""

from __future__ import annotations

from ns_exploration.conjectures.l0046_greedy_fullsym import SUPPORT_CR0011
from ns_exploration.conjectures.l0047_greedy_fullsym import (
    SUPPORT_CR0012,
    cr0012_record,
    lemma_l0047,
)
from ns_exploration.validation.l0047_certificate import (
    build_l0047_certificate,
    verify_l0047_certificate,
)


def test_l0047_cr0012_and_certificate():
    assert set(SUPPORT_CR0011).issubset(set(SUPPORT_CR0012))
    assert len(SUPPORT_CR0012) == 41
    assert 147 in SUPPORT_CR0012
    b = lemma_l0047(n=24)
    assert b.support_cr0012 == list(SUPPORT_CR0012)
    assert b.C_fullsym_cr0012 <= b.C_dagger
    assert b.closes_cr0012
    assert not b.closes_c0007_all_ic
    cr = cr0012_record(b)
    assert cr["id"] == "C-R-0012"
    cert = build_l0047_certificate(n=24, bound=b)
    ok, checks = verify_l0047_certificate(cert.as_dict())
    assert ok and all(checks.values())
