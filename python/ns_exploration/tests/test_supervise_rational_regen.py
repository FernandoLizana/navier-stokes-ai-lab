"""Tests for external regen supervisor helpers."""

from __future__ import annotations

import time
from pathlib import Path

from scripts.supervise_c0008_rational_regen import _heartbeat_age_sec, _paths, _tag


def test_tag_rational():
    assert _tag(24, True) == "n24_rational"
    assert _tag(24, False) == "n24"


def test_paths_contain_stop_and_heartbeat():
    p = _paths(24, True)
    assert "supervisor_stop" in p["stop_flag"].name
    assert "heartbeat" in p["heartbeat"].name


def test_heartbeat_age(tmp_path: Path):
    hb = tmp_path / "hb.txt"
    assert _heartbeat_age_sec(hb) is None
    hb.write_text("x", encoding="utf-8")
    age = _heartbeat_age_sec(hb)
    assert age is not None
    assert 0 <= age < 2
    time.sleep(0.05)
    assert _heartbeat_age_sec(hb) > age
