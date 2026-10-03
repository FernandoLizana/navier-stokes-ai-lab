"""Tests for L-0018 spectral enstrophy envelope / C-S-0002 closure."""

from __future__ import annotations

from ns_exploration.conjectures.l0018_envelope import (
    dealias_K2,
    envelope_Omega,
    lemma_l0018,
    max_K2_upto,
)
from ns_exploration.validation.l0018_certificate import (
    build_l0018_envelope_certificate,
    verify_l0018_envelope_certificate,
)


def test_l0018_closes_cs0002():
    b = lemma_l0018()
    assert b.closes_cs0002
    assert b.K2_worst == 75
    assert b.n_worst == 16
    assert b.Omega_cap == 37.5
    assert b.Omega_cap < b.cs0002_M


def test_l0018_K2_table():
    assert dealias_K2(16) == 75
    assert dealias_K2(12) == 27
    K2, n = max_K2_upto(16)
    assert (K2, n) == (75, 16)
    assert envelope_Omega(16, 0.5) == 37.5


def test_l0018_certificate():
    cert = build_l0018_envelope_certificate()
    assert cert.closes
    assert cert.Omega_cap_hi <= cert.R
    ok, checks = verify_l0018_envelope_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_l0018_closes_c0003_equality():
    """All-IC C-0003 takes M = K² E0 = 37.5 exactly."""
    cert = build_l0018_envelope_certificate(R=37.5)
    assert cert.closes
    assert abs(cert.Omega_cap_hi - 37.5) < 1e-12


def test_l0018_does_not_close_c0002_scale():
    """C-0002 had M≈20.32; envelope 37.5 is too weak (and C-0002 is refuted)."""
    b = lemma_l0018()
    assert b.Omega_cap > 20.32


def test_l0018_tamper_fails():
    cert = build_l0018_envelope_certificate().as_dict()
    cert["Omega_cap_hi"] = cert["R"] + 1.0
    ok, checks = verify_l0018_envelope_certificate(cert)
    assert not ok
    assert not checks["Omega_le_R"]
