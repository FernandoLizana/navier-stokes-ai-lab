"""Tests: Stokes max-mode refutes C-0002; C-0003 closed by L-0018."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0002_stokes_refuter import (
    stokes_Omega_exact,
    stokes_max_mode_ic,
)
from ns_exploration.conjectures.l0018_envelope import lemma_l0018
from ns_exploration.diagnostics.metrics import compute_diagnostics
from ns_exploration.spectral.integrators import step_etd_rk2, step_rk4
from ns_exploration.validation.l0018_certificate import build_l0018_envelope_certificate


def test_stokes_refuter_exceeds_c0002_M():
    M = 20.320076249551253
    u0, meta = stokes_max_mode_ic(16, 0.5)
    assert meta["K2"] == 75
    assert abs(compute_diagnostics(u0, 0.1).enstrophy - 37.5) < 1e-10
    Om_ex = stokes_Omega_exact(75, 0.5, 0.1, 0.02)
    assert Om_ex > M
    assert abs(Om_ex - 27.78068327556442) < 1e-8


def test_stokes_etd_matches_exact():
    u0, meta = stokes_max_mode_ic(16, 0.5)
    u = u0.copy()
    for _ in range(20):
        u = step_etd_rk2(u, 0.1, 1e-3)
    Om = compute_diagnostics(u, 0.1).enstrophy
    Om_ex = stokes_Omega_exact(meta["K2"], 0.5, 0.1, 0.02)
    assert abs(Om - Om_ex) / Om_ex < 1e-10
    u_rk = u0.copy()
    for _ in range(20):
        u_rk = step_rk4(u_rk, 0.1, 1e-3)
    Om_rk = compute_diagnostics(u_rk, 0.1).enstrophy
    assert abs(Om - Om_rk) / Om < 1e-8


def test_c0003_proved_by_envelope():
    b = lemma_l0018()
    assert b.Omega_cap == 37.5
    cert = build_l0018_envelope_certificate(R=37.5)
    assert cert.closes
