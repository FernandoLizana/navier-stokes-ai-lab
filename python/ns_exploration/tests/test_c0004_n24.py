"""Tests for C-0004 / C-0005 (both proved on N≤24)."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.c0004_n24_line import (
    build_c0004_and_c0005,
    stokes_floor_N,
)
from ns_exploration.conjectures.l0018_envelope import lemma_l0018
from ns_exploration.validation.l0018_certificate import build_l0018_envelope_certificate


def test_stokes_floor_n24():
    floor, K2 = stokes_floor_N(24, 0.5, 0.1, 0.02)
    assert K2 == 147
    assert abs(floor - 40.82462307930434) < 1e-8


def test_c0004_closed_by_envelope():
    b = lemma_l0018(n_max=24, cs0002_M=73.5)
    assert b.K2_worst == 147
    assert b.Omega_cap == 73.5
    assert b.closes_cs0002
    cert = build_l0018_envelope_certificate(n_max=24, R=73.5)
    assert cert.closes


def test_c0005_between_stokes_and_envelope_and_proved():
    out = build_c0004_and_c0005()
    assert out["stokes_floor"] < out["C0005_M"] < out["C0004_M"]
    assert out["C0004_proved"]
    c5 = json.loads(
        Path("conjectures/proved_restricted/C-0005.json").read_text(encoding="utf-8")
    )
    c4 = json.loads(
        Path("conjectures/proved_restricted/C-0004.json").read_text(encoding="utf-8")
    )
    assert c5["id"] == "C-0005" and c5["status"] == "proved_restricted"
    assert c5["lemma"] == "L-0024"
    assert c4["id"] == "C-0004" and c4["status"] == "proved_restricted"
    assert not Path("conjectures/active/C-0005.json").exists()


def test_l0019_stokes_majorant():
    from ns_exploration.conjectures.l0019_stokes_majorant import lemma_l0019
    from ns_exploration.validation.l0019_certificate import (
        build_l0019_certificate,
        verify_l0019_certificate,
    )

    b = lemma_l0019()
    assert b.K2_star == 147
    assert abs(b.Stokes_floor - 40.82462307930434) < 1e-8
    assert b.below_c0005
    cert = build_l0019_certificate()
    ok, checks = verify_l0019_certificate(cert.as_dict())
    assert ok and all(checks.values())
