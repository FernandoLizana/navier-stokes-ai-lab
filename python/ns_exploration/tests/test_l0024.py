"""Tests for L-0024 spectral-defect majorant and C-0005 closure."""

from __future__ import annotations

from pathlib import Path

from ns_exploration.conjectures.l0024_spectral_defect import (
    dealias_shell_stats,
    defect_ode_omega_T,
    lemma_l0024,
    max_shell_has_no_self_triads,
    max_shell_modes,
)
from ns_exploration.validation.l0024_certificate import (
    build_l0024_certificate,
    verify_l0024_certificate,
)


def test_max_shell_geometry():
    modes = max_shell_modes(24)
    assert len(modes) == 8
    assert all(a * a + b * b + c * c == 147 for a, b, c in modes)
    assert max_shell_has_no_self_triads(24)
    s = dealias_shell_stats(24)
    assert s["K2"] == 147 and s["gap"] == 13 and s["m_K"] == 8


def test_stokes_endpoint_matches_floor():
    s = dealias_shell_stats(24)
    om = defect_ode_omega_T(0.5, 73.5, 0.1, 0.02, s)
    assert abs(om - 40.82462307930434) < 1e-6


def test_l0024_closes_c0005():
    b = lemma_l0024(n_grid=41)
    assert b.closes_c0005
    assert b.Omega_T_worst < b.c0005_M
    assert b.Omega_T_worst < 45.0


def test_l0024_certificate():
    cert = build_l0024_certificate()
    ok, checks = verify_l0024_certificate(cert.as_dict())
    assert ok and all(checks.values())
