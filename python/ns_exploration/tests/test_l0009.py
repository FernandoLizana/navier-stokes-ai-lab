"""Tests for L-0009 two-scale bound and N5 certificate."""

from __future__ import annotations

from ns_exploration.conjectures.l0009_twoscale import lemma_l0009, shell_mode_count
from ns_exploration.validation.l0009_certificate import (
    build_l0009_certificate,
    verify_l0009_certificate,
)


def test_shell_mode_count_eucl_k2():
    assert shell_mode_count(12, 2, "euclidean") == 33


def test_l0009_beats_l0008_at_n12():
    b = lemma_l0009(n=12, k_ic=2, ic_kind="euclidean", cs0002_M=None)
    assert b.two_scale_closes
    assert b.which_best == "L-0009-two-scale"
    assert b.two_scale_cap is not None
    assert b.two_scale_cap < b.L0008_cap
    assert 100.0 < b.two_scale_cap < 110.0


def test_l0009_cs0002_class_does_not_close():
    b = lemma_l0009(n=16, k_ic=4, ic_kind="linf")
    assert b.two_scale_closes is False
    assert b.proves_cs0002 is False
    assert b.best_cap > 1000.0


def test_l0009_certificate():
    cert = build_l0009_certificate()
    ok, checks = verify_l0009_certificate(cert.as_dict())
    assert ok and all(checks.values())
    b = lemma_l0009(n=12, k_ic=2, ic_kind="euclidean", cs0002_M=None)
    assert cert.Omega_bound_hi >= float(b.two_scale_cap) - 1e-2


def test_l0009_rejects_understated():
    cert = build_l0009_certificate().as_dict()
    cert["Omega_bound_hi"] = 1.0
    ok, checks = verify_l0009_certificate(cert)
    assert not ok
    assert not checks["bound_is_valid_upper"]
