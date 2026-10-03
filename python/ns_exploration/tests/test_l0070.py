"""Tests for L-0070 all-IC C-0007 closure (Option A)."""

from __future__ import annotations

from pathlib import Path

from ns_exploration.conjectures.l0070_c0007_all_ic_l0024 import (
    build_c0007_proved,
    build_c0008_sharp,
    lemma_l0070,
)
from ns_exploration.validation.l0070_certificate import (
    build_l0070_certificate,
    verify_l0070_certificate,
)


def test_l0070_closes_all_ic_at_l0024_majorant():
    b = lemma_l0070(n_grid=41)
    assert b.closes_c0007_all_ic
    assert not b.closes_c0008_sharp
    assert b.monotone_in_Omega0
    assert abs(b.proved_bound_M - 41.74337593971423) < 1e-6
    assert abs(b.sharp_target_M - 41.283999509509286) < 1e-6
    assert b.Omega0_worst == b.E0


def test_l0070_certificate():
    cert = build_l0070_certificate()
    ok, checks = verify_l0070_certificate(cert.as_dict())
    assert ok and all(checks.values())


def test_c0007_proved_json_shape():
    b = lemma_l0070()
    c7 = build_c0007_proved(b)
    assert c7["status"] == "proved_restricted"
    assert c7["lemma"] == "L-0070"
    assert c7["sharp_refinement"] == "C-0008"
    assert c7["proved_bound_M"] < 42.0


def test_c0008_sharp_open():
    b = lemma_l0070()
    c8 = build_c0008_sharp(b)
    assert c8["status"] == "exploring"
    assert c8["parent"] == "C-0007"
    assert c8["parent_proved_M"] == b.proved_bound_M
