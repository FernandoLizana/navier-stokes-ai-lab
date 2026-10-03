"""Tests for L-0042 Sym-strengthening techo."""

from __future__ import annotations

from ns_exploration.conjectures.l0042_sym_strengthen_techo import lemma_l0042
from ns_exploration.validation.l0042_certificate import (
    build_l0042_certificate,
    verify_l0042_certificate,
)


def test_l0042_strengthenings_fail_123456():
    b = lemma_l0042(n=24)
    assert b.C_sym_123456 > b.C_dagger
    assert b.C_triangle_123456 > b.C_sym_123456
    assert b.C_hybrid_AB > b.C_dagger
    assert b.C_rank1_lo_123456 < b.C_dagger
    assert not b.closes_123456
    assert not b.closes_c0007_all_ic


def test_l0042_certificate():
    cert = build_l0042_certificate(n=24)
    ok, checks = verify_l0042_certificate(cert.as_dict())
    assert ok and all(checks.values())
