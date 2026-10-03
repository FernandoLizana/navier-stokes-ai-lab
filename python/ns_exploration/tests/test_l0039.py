"""Tests for L-0039 sparse support Shor and C-R-0005."""

from __future__ import annotations

from ns_exploration.conjectures.l0038_band_svd_shor import exact_band_C_shor
from ns_exploration.conjectures.l0039_sparse_support_shor import (
    QUAD_1346,
    SUPPORT_A,
    build_cr0005,
    lemma_l0039,
    sparse_band_C_shor,
)
from ns_exploration.validation.l0039_certificate import (
    build_l0039_certificate,
    verify_l0039_certificate,
)


def test_l0039_quad_and_supports():
    b = lemma_l0039(n=24)
    assert b.quad_1346_C <= b.C_dagger
    assert b.support_A_C <= b.C_dagger
    assert b.support_B_C <= b.C_dagger
    assert len(b.support_A or []) == 14
    assert b.closes_cr0005
    assert not b.closes_c0007_all_ic


def test_l0039_sparse_matches_dense_quad():
    dense = exact_band_C_shor(24, QUAD_1346)["C_shor"]
    sparse = sparse_band_C_shor(24, QUAD_1346)["C_shor"]
    assert abs(dense - sparse) < 1e-8


def test_cr0005_build():
    cr = build_cr0005(n=24)
    assert cr["id"] == "C-R-0005"
    assert cr["shells"] == list(SUPPORT_A)
    assert cr["C_shor"] <= cr["C_dagger"]


def test_l0039_certificate():
    cert = build_l0039_certificate(n=24)
    ok, checks = verify_l0039_certificate(cert.as_dict())
    assert ok and all(checks.values())
