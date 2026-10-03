"""Tests for L-0036 signed Shor stretch techo."""

from __future__ import annotations

from ns_exploration.conjectures.l0036_signed_shor_stretch import lemma_l0036
from ns_exploration.validation.l0036_certificate import (
    build_l0036_certificate,
    verify_l0036_certificate,
)


def test_l0036_shor_exceeds_cdagger():
    b = lemma_l0036(n=12, iters=10, include_n24_probe=False)
    assert b.C_shor_lower > b.C_dagger
    assert b.L_op_lower > 0.0
    assert not b.closes_c0007


def test_l0036_certificate():
    cert = build_l0036_certificate(n=12, iters=10)
    ok, checks = verify_l0036_certificate(cert.as_dict())
    assert ok and all(checks.values())
