"""Tests for L-0021 mono-radial quartic Gram bound."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.conjectures.l0021_quartic_shell import lemma_l0021
from ns_exploration.validation.l0021_certificate import (
    build_l0021_certificate,
    verify_l0021_certificate,
)


def test_l0021_beats_envelope_not_c0005():
    b = lemma_l0021(iters=40)
    assert b.alpha_star == 14.0 or abs(b.alpha_star - 14.0) < 1e-6
    assert b.r_star == 74
    assert b.beats_envelope
    assert not b.closes_c0005
    assert b.Omega_radial_H1 < b.envelope_cap
    assert b.Omega_radial_H1 > b.c0005_M


def test_l0021_certificate_and_cr0001():
    cert = build_l0021_certificate(iters=40)
    ok, checks = verify_l0021_certificate(cert.as_dict(), iters=40)
    assert ok and all(checks.values())
    # sprint artifact may exist after run; if not, still ok
    p = Path("conjectures/proved_restricted/C-R-0001.json")
    if p.exists():
        cr = json.loads(p.read_text(encoding="utf-8"))
        assert cr["id"] == "C-R-0001"
        assert cr["status"] == "proved_restricted"
        assert cr["proposed_bound_M"] < 73.5
