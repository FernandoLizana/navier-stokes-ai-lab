"""Tests for L-0037 shell-diff / bi-shell Shor techo."""

from __future__ import annotations

from ns_exploration.conjectures.l0037_shell_diff_stretch import (
    bishell_shor_C_lower,
    lemma_l0037,
    mono_radial_stretch_identity,
)
from ns_exploration.validation.l0037_certificate import (
    build_l0037_certificate,
    verify_l0037_certificate,
)


def test_l0037_mono_identity_documented():
    assert "stretch" in mono_radial_stretch_identity().lower()
    assert "0" in mono_radial_stretch_identity()


def test_l0037_near_meets_far_fails():
    c_near = bishell_shor_C_lower(24, 1, 2, iters=14)["C_shor_lower"]
    c_far = bishell_shor_C_lower(24, 2, 134, iters=14)["C_shor_lower"]
    b = lemma_l0037(n=24, iters=12, max_pairs=25)
    assert c_near <= b.C_dagger * 1.05
    assert c_far > b.C_dagger
    assert b.C_bishell_shor_max_sample > b.C_dagger
    assert not b.closes_c0007


def test_l0037_certificate():
    cert = build_l0037_certificate(n=24)
    ok, checks = verify_l0037_certificate(cert.as_dict())
    assert ok and all(checks.values())
