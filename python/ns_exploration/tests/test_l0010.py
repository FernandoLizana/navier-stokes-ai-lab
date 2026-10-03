"""Tests for L-0010 short-time Duhamel bound and N5 certificate."""

from __future__ import annotations

from ns_exploration.conjectures.l0010_duhamel_short import (
    duhamel_bootstrap_closes,
    lemma_l0010,
)
from ns_exploration.validation.l0010_certificate import (
    build_l0010_short_certificate,
    verify_l0010_short_certificate,
)


def test_duhamel_closes_short_not_long():
    b = lemma_l0010(n=16, k_ic=4, ic_kind="linf")
    assert b.T_duhamel_for_cs0002 > 0.0005
    assert b.T_duhamel_for_cs0002 < 0.005
    assert b.duhamel_closes_at_T_target is False
    assert b.proves_cs0002_on_short_horizon is True


def test_bootstrap_at_Tstar():
    b = lemma_l0010()
    K = b.K_full_squared**0.5
    ok, eps, eh = duhamel_bootstrap_closes(
        b.B, K, b.U_L, b.W, b.cs0002_M, b.T_duhamel_for_cs0002
    )
    assert ok
    assert eh <= eps + 1e-12


def test_l0010_certificate():
    cert = build_l0010_short_certificate()
    assert cert.bootstrap_holds
    assert cert.T_cert < 0.02
    assert cert.T_cert < cert.T_star
    ok, checks = verify_l0010_short_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0010_rejects_bad_T():
    cert = build_l0010_short_certificate().as_dict()
    cert["T_cert"] = 0.02  # too long
    ok, checks = verify_l0010_short_certificate(cert)
    assert not ok
    assert not checks["bootstrap_holds"]
