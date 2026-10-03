"""Tests for L-0007 (shell IC → full-grid) and its N5 exponential certificate."""

from __future__ import annotations

from ns_exploration.conjectures.l0007_shell_ic_fullgrid import (
    lemma_l0007,
    linf_shell_K2,
    euclidean_shell_K2,
)
from ns_exploration.validation.l0007_certificate import (
    build_l0007_exp_certificate,
    verify_l0007_exp_certificate,
)


def test_linf_shell_K2():
    assert linf_shell_K2(4) == 48  # 3*16
    assert linf_shell_K2(2) == 12


def test_euclidean_shell_K2_n12():
    assert euclidean_shell_K2(12, 2) == 4
    assert euclidean_shell_K2(12, 3) == 9


def test_l0007_improves_on_full_l0003_at_n12():
    """Shell Ω0 makes the exponential beat the uniform full-grid ceiling at N=12."""
    b = lemma_l0007(n=12, k_ic=2, ic_kind="euclidean", cs0002_M=None)
    assert b.Omega0_cap == 2.0  # K²E0 = 4*0.5
    assert b.best_cap < b.L0003_cap
    assert b.which_best == "L-0007-exp"
    assert not b.proves_cs0002  # still far from conjectured ~48


def test_l0007_cs0002_class_does_not_close_to_M():
    b = lemma_l0007(n=16, k_ic=4, ic_kind="linf")
    assert b.K_ic_squared == 48
    assert b.M_full == 1331
    assert b.best_cap > float(b.cs0002_M)
    assert b.proves_cs0002 is False


def test_l0007_exp_certificate_verifies():
    cert = build_l0007_exp_certificate(n=12, k_ic=2, ic_kind="euclidean")
    ok, checks = verify_l0007_exp_certificate(cert.as_dict())
    assert ok
    assert all(checks.values())
    # Must enclose the lemma's float exponential
    b = lemma_l0007(n=12, k_ic=2, ic_kind="euclidean", cs0002_M=None)
    assert cert.Omega_bound_hi >= b.L0001_cap - 1e-6


def test_l0007_exp_rejects_understated():
    cert = build_l0007_exp_certificate().as_dict()
    cert["Omega_bound_hi"] = 1.0
    ok, checks = verify_l0007_exp_certificate(cert)
    assert not ok
    assert not checks["bound_is_valid_upper"]
