"""L-0075 hybrid low-slab sprint tests."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.skipif(
    not (ROOT / "experiments/terminal_weighted/hybrid_low_slab_summary.json").is_file(),
    reason="hybrid summary missing",
)
def test_hybrid_low_slab_summary():
    s = json.loads(
        (ROOT / "experiments/terminal_weighted/hybrid_low_slab_summary.json").read_text(
            encoding="utf-8"
        )
    )
    assert s["closes_c0008_sharp"] is False
    assert s["best_all_ic_low_slab_omega_T"] > s["M_sharp_c0008"]
    assert s["best_label"] in s["candidates"]
    assert s["best_low_slab_omega_T"] == min(s["candidates"].values())


@pytest.mark.skipif(
    not (ROOT / "conjectures/active/L-0075.json").is_file(),
    reason="L-0075 registry missing",
)
def test_l0075_registry():
    reg = json.loads((ROOT / "conjectures/active/L-0075.json").read_text(encoding="utf-8"))
    assert reg["lemma_id"] == "L-0075"
    assert reg["closes_c0008_sharp"] is False
