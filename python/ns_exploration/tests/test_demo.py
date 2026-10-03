"""Tests for the public demo entrypoint."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from ns_exploration.demo import build_parser, main, run_demo


def test_demo_rejects_odd_n():
    with pytest.raises(SystemExit):
        args = build_parser().parse_args(["--n", "7"])
        run_demo(args)


def test_demo_rejects_negative_nu():
    with pytest.raises(SystemExit):
        args = build_parser().parse_args(["--nu", "-0.1"])
        run_demo(args)


def test_demo_rejects_too_many_steps():
    with pytest.raises(SystemExit):
        args = build_parser().parse_args(
            ["--n", "8", "--dt", "1e-6", "--t-end", "1.0"]
        )
        run_demo(args)


def test_demo_completes_and_writes_outputs(tmp_path: Path):
    out = tmp_path / "run"
    args = build_parser().parse_args(
        [
            "--n",
            "8",
            "--nu",
            "0.05",
            "--dt",
            "0.001",
            "--t-end",
            "0.005",
            "--out-dir",
            str(out),
        ]
    )
    summary = run_demo(args)
    assert summary["result"]["finished_normally"] is True
    assert summary["result"]["steps"] == 5
    assert (out / "diagnostics.csv").is_file()
    assert (out / "config.json").is_file()
    assert (out / "summary.json").is_file()
    assert (out / "diagnostics.png").is_file()
    data = json.loads((out / "summary.json").read_text(encoding="utf-8"))
    assert data["evidence_level"] == "N2"
    assert data["result"]["energy_final"] > 0
    assert data["result"]["div_l2_final"] < 1e-10


def test_demo_main_exit_code(tmp_path: Path):
    code = main(
        [
            "--n",
            "8",
            "--dt",
            "0.001",
            "--t-end",
            "0.005",
            "--out-dir",
            str(tmp_path / "cli"),
        ]
    )
    assert code == 0
