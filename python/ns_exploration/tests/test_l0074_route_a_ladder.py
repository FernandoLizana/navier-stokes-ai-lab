"""L-0074 Route A ladder negative result tests."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.skipif(
    not (ROOT / "experiments/terminal_weighted/route_a_ladder_full.json").is_file(),
    reason="ladder summary missing",
)
def test_l0074_ladder_equal_12_optimal():
    from ns_exploration.validation.l0074_certificate import verify_l0074_ladder

    ok, msgs = verify_l0074_ladder()
    assert ok, msgs


@pytest.mark.skipif(
    not (ROOT / "conjectures/active/L-0074.json").is_file(),
    reason="L-0074 registry missing",
)
def test_l0074_registry():
    reg = json.loads((ROOT / "conjectures/active/L-0074.json").read_text(encoding="utf-8"))
    assert reg["status"] == "partition_ladder_refuted"
    assert reg["best_label"] == "equal_12"
    assert reg["closes_c0008"] is False
    assert reg["best_omega_T_hi"] > reg["M_target"]
