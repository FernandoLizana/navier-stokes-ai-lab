"""Tests for C-0008 repair finalize pipeline."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.validation.c0008_repair_pipeline import run_finalize_pipeline

REPO = Path(__file__).resolve().parents[3]
BAND_MANIFEST = REPO / "experiments/terminal_weighted/shell_manifest_repair_band6_one_inf.json"


def test_finalize_band_pipeline_verifies():
    if not BAND_MANIFEST.is_file():
        return
    rep = run_finalize_pipeline(
        BAND_MANIFEST,
        prec=128,
        workers=1,
        out_cert=REPO / "certificates/CERT-C0008-repair-band6-one-inf.json",
        out_cluster=REPO / "experiments/terminal_weighted/cluster_manifest_repair_band6_one_inf.json",
    )
    assert rep["verify_ok"], rep.get("issues")
    assert rep["best_route"] in ("A", "A_prime", "B")
    assert "cluster_manifest" in rep
