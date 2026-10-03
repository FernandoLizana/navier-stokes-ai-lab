"""Tests for L-0025 / C-0006 / C-0007."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0024_spectral_defect import max_shell_has_no_self_triads
from ns_exploration.conjectures.l0025_multi_n_defect import c0007_open_target, lemma_l0025
from ns_exploration.validation.l0025_certificate import (
    build_l0025_certificate,
    verify_l0025_certificate,
)


def test_l0025_n32_closes_c0006():
    assert max_shell_has_no_self_triads(32)
    b = lemma_l0025(n=32, n_grid=41)
    assert b.K2 == 300
    assert b.closes_c0006
    assert b.Omega_T_worst < b.c0006_M
    assert b.Stokes_floor < b.Omega_T_worst


def test_l0025_certificate():
    cert = build_l0025_certificate(n=32)
    ok, checks = verify_l0025_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_c0007_gap_is_nontrivial():
    g = c0007_open_target(24)
    assert g["Stokes_floor"] < g["proposed_bound_M"] < g["L0024_majorant"]
    assert g["open_gap"] > 0.5
