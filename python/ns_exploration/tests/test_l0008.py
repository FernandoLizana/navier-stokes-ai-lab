"""Tests for L-0008 cascade bound and N5 cosh certificate."""

from __future__ import annotations

from ns_exploration.conjectures.l0008_cascade import cascade_cosh_bound, lemma_l0008
from ns_exploration.validation.l0008_certificate import (
    build_l0008_cosh_certificate,
    verify_l0008_cosh_certificate,
)


def test_cascade_improves_on_l0007_at_n12():
    b = lemma_l0008(n=12, k_ic=2, ic_kind="euclidean", cs0002_M=None)
    assert b.which_best == "L-0008-cosh"
    assert b.cosh_cap < b.L0007_cap
    assert b.cosh_cap < 500.0  # ~394


def test_l0008_does_not_prove_cs0002_at_n16():
    b = lemma_l0008(n=16, k_ic=4, ic_kind="linf")
    assert b.proves_cs0002 is False
    # cosh explodes; best falls back to L-0003/L-0007 ~2e5
    assert b.best_cap > 1000.0


def test_cosh_formula_at_t0():
    A, alpha, cap = cascade_cosh_bound(Omega0_shell=2.0, M_full=343, K_full=5.196, E0=0.5, t=0.0)
    assert abs(cap - 2.0) < 1e-9


def test_l0008_certificate_verifies():
    cert = build_l0008_cosh_certificate(n=12, k_ic=2, ic_kind="euclidean")
    ok, checks = verify_l0008_cosh_certificate(cert.as_dict())
    assert ok and all(checks.values())
    b = lemma_l0008(n=12, k_ic=2, ic_kind="euclidean", cs0002_M=None)
    assert cert.Omega_bound_hi >= b.cosh_cap - 1e-3


def test_l0008_rejects_understated():
    cert = build_l0008_cosh_certificate().as_dict()
    cert["Omega_bound_hi"] = 1.0
    ok, checks = verify_l0008_cosh_certificate(cert)
    assert not ok
    assert not checks["bound_is_valid_upper"]
