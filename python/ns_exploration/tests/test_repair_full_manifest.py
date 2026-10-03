"""Tests for full-dealias repair manifest pipeline."""

from __future__ import annotations

import json

import pytest
from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.repair_full_manifest import (
    _is_memory_error,
    assemble_repair_manifest,
    build_repair_full_one_inf_manifest,
    pick_workers,
)
from ns_exploration.validation.c0008_repair_band_cert import build_repair_full_partial_certificate
from ns_exploration.validation.c0008_verify import verify_terminal_certificate


def test_assemble_repair_manifest_structure():
    fake = {
        1: {
            "shell": 1,
            "n": 12,
            "D": 52,
            "C_term_hi": "1.0",
            "method": "one_inf_mpfr_up",
            "n_updates": 10,
            "evidence_level": "N4",
            "payload_sha256": "abc",
        }
    }
    m = assemble_repair_manifest(
        n=12, radii=(1, 2), blocks_by_shell=fake, prec=128, complete=False
    )
    assert m["bound_method"] == "one_inf_only"
    assert m["repair_complete"] is False
    assert len(m["blocks"]) == 1
    assert "manifest_sha256" in m


def test_repair_full_partial_certificate_from_manifest():
    fake = assemble_repair_manifest(
        n=12,
        radii=(1, 2),
        blocks_by_shell={
            1: {
                "shell": 1,
                "n": 12,
                "D": 52,
                "C_term_hi": "3.5",
                "method": "one_inf_mpfr_up",
                "n_updates": 10,
                "evidence_level": "N4",
                "payload_sha256": "x",
                "C_term_one_inf_hi": "3.5",
            }
        },
        prec=128,
        complete=False,
    )
    cert = build_repair_full_partial_certificate(fake, prec=128)
    ok, issues = verify_terminal_certificate(
        cert, manifest=fake, manifest_path=None, prec=128
    )
    assert ok, issues
    assert cert["scope"] == "repair_full_partial"


def test_build_repair_full_n12_shell1():
    radii = tuple(all_dealias_radii(12))
    m = build_repair_full_one_inf_manifest(
        12, prec=128, workers=1, resume=True, shell_limit=1
    )
    assert m["n"] == 12
    assert len(m["blocks"]) >= 1
    assert m["blocks"][0]["method"] == "one_inf_mpfr_up"
    assert float(mpfr(m["blocks"][0]["C_term_hi"])) > 0


def test_pick_workers_respects_caps():
    w = pick_workers(99, n_pending=3, ram_gb_per_worker=100.0)
    assert 1 <= w <= 3


def test_is_memory_error():
    assert _is_memory_error(MemoryError("x"))
    assert _is_memory_error(RuntimeError("Out of memory allocating tensor"))
    assert not _is_memory_error(ValueError("bad shell"))
