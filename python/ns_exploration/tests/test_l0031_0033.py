"""Tests for L-0031/32/33 techos and C-R-0003 seeds."""

from __future__ import annotations

from ns_exploration.conjectures.l0031_0033_techos import lemma_l0031_33
from ns_exploration.validation.l0031_0033_certificate import (
    build_l0031_33_certificate,
    verify_l0031_33_certificate,
)


def test_l0031_33_techos_and_zero_alpha():
    b = lemma_l0031_33(n=24, n_grid=15)
    assert b.compatible_worst_OmT > b.c0007_M
    assert b.compatible_worst_OmT >= b.L0024_majorant - 0.05
    assert not b.Smax_beats_young
    assert b.C_from_cubic_F > b.C_dagger
    assert b.zero_alpha_shells and 147 in b.zero_alpha_shells
    assert b.zero_alpha_duhamel_OmT <= b.c0007_M + 1e-9
    assert not b.closes_c0007


def test_l0031_33_certificate():
    cert = build_l0031_33_certificate(n=24)
    ok, checks = verify_l0031_33_certificate(cert.as_dict())
    assert ok and all(checks.values())
