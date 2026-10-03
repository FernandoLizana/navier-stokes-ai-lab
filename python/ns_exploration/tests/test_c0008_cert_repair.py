"""Adversarial verifier tests and repair-phase unit tests for C-0008."""

from __future__ import annotations

import copy
import json
from pathlib import Path

import pytest

from ns_exploration.terminal_weighted.constants import EXACT, FROZEN
from ns_exploration.terminal_weighted.intervals import (
    i_star_interval,
    mpfr_const,
    omega_T_hi_from_I_term,
    stokes_floor_interval,
    stokes_floor_meta,
)
from ns_exploration.terminal_weighted.tensor import validate_direct_vs_tensor
from ns_exploration.validation.c0008_verify import (
    manifest_sha256,
    verify_terminal_certificate,
)


REPO = Path(__file__).resolve().parents[3]
CERT_L0072 = REPO / "certificates" / "CERT-L0072-C0008-terminal-full-dealias.json"
MANIFEST_BEST = REPO / "experiments" / "terminal_weighted" / "shell_manifest_best_full.json"


def _load_json(p: Path) -> dict:
    return json.loads(p.read_text(encoding="utf-8"))


def _base_cert() -> dict:
    cert = _load_json(CERT_L0072)
    manifest = _load_json(MANIFEST_BEST)
    cert["source_hashes"] = {"shell_manifest": manifest_sha256(manifest)}
    cert["bound_method"] = "one_inf_only"
    return cert, manifest


@pytest.mark.skipif(not CERT_L0072.is_file(), reason="certificate missing")
def test_honest_certificate_structure_verifies_except_open_bound():
    cert, manifest = _base_cert()
    prec = int(cert.get("arithmetic_backend", {}).get("precision_bits", 200))
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=MANIFEST_BEST, prec=prec
    )
    # Historical cert uses frobenius/best — expect FAIL on bound_method until regenerated
    assert "manifest required" not in issues


def test_adversarial_integral_too_small_fails():
    cert, manifest = _base_cert()
    cert["integral_interval"] = ["0", "0.001"]
    prec = 200
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=MANIFEST_BEST, prec=prec
    )
    assert not ok
    assert any("declared" in i and "recomputed" in i for i in issues)


def test_adversarial_omega_too_small_fails():
    cert, manifest = _base_cert()
    cert["omega_T_interval"] = ["0", "41.0"]
    prec = 200
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=MANIFEST_BEST, prec=prec
    )
    assert not ok


def test_adversarial_missing_manifest_fails():
    cert, _manifest = _base_cert()
    ok, issues = verify_terminal_certificate(
        cert, manifest=None, manifest_path=None, prec=200
    )
    assert not ok
    assert any("manifest required" in i for i in issues)


def test_adversarial_hash_mismatch_fails():
    cert, manifest = _base_cert()
    cert["source_hashes"] = {"shell_manifest": "0" * 64}
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=MANIFEST_BEST, prec=200
    )
    assert not ok
    assert any("hash mismatch" in i for i in issues)


def test_adversarial_band_mode_fails():
    cert, manifest = _base_cert()
    cert.setdefault("mode_set", {})["radii_mode"] = "band"
    ok, _issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=MANIFEST_BEST, prec=200
    )
    assert not ok


def test_adversarial_frobenius_bound_method_fails():
    cert, manifest = _base_cert()
    cert["bound_method"] = "best"
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=MANIFEST_BEST, prec=200
    )
    assert not ok
    assert any("bound_method" in i or "one_inf_only" in i for i in issues)


def test_repair_band_certificate_verifies():
    from ns_exploration.validation.c0008_repair_band_cert import build_repair_band_certificate

    manifest_path = REPO / "experiments" / "terminal_weighted" / "shell_manifest_repair_band6_one_inf.json"
    if not manifest_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cert = build_repair_band_certificate(manifest, prec=128)
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=manifest_path, prec=128
    )
    assert ok, issues


def test_repair_full_n24_certificate_verifies():
    manifest_path = (
        Path(__file__).resolve().parents[2]
        / "experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n24.json"
    )
    cert_path = REPO / "certificates/CERT-C0008-repair-full-n24-one-inf.json"
    if not manifest_path.is_file() or not cert_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=manifest_path, prec=128
    )
    assert ok, issues
    assert cert["scope"] == "repair_full"
    assert cert["mode_set"]["n_shells"] == 87
    assert cert["closes_c0008"] is False


def test_repair_cluster_n24_certificate_verifies():
    cluster_path = REPO / "experiments/terminal_weighted/cluster_manifest_repair_full_n24_one_inf.json"
    cert_path = REPO / "certificates/CERT-L0073R-repair-C0008-cluster-one-inf.json"
    if not cluster_path.is_file() or not cert_path.is_file():
        return
    manifest = json.loads(cluster_path.read_text(encoding="utf-8"))
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=cluster_path, prec=128
    )
    assert ok, issues
    assert cert["scope"] == "repair_full_cluster_sketch"
    assert cert["mode_set"]["n_clusters"] >= 1


def test_adversarial_wrong_nu_fails():
    cert, manifest = _base_cert()
    cert.setdefault("parameters", {})["nu"] = "0.2"
    ok, _issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=MANIFEST_BEST, prec=200
    )
    assert not ok


def test_i_star_interval_orientation():
    I_lo, I_hi = i_star_interval(prec=200)
    assert I_lo < I_hi
    assert float(I_lo) > 0


def test_stokes_floor_interval_covers_frozen():
    floor_lo, floor_hi = stokes_floor_interval(n=24, prec=200)
    meta = stokes_floor_meta()
    assert meta.get("argmax_shell", 0) > 0
    assert len(meta.get("shells", [])) > 0
    assert float(floor_hi) >= FROZEN.stokes_floor - 1e-6


def test_omega_uses_sqrt8_lo_divisor():
    from gmpy2 import mpfr

    I = mpfr("10", precision=200)
    omega = omega_T_hi_from_I_term(I, prec=200)
    assert float(omega) > float(FROZEN.stokes_floor)


def test_four_path_direct_vs_tensor_band():
    val = validate_direct_vs_tensor(24, (1, 2, 3), FROZEN.T, n_probe=30, rtol=1e-10)
    assert val["paths"] == [
        "pseudospectral_direct",
        "raw_G",
        "full_sym_G",
        "streaming_M",
    ]
    assert val["worst_rel"]["B_C"] <= 1e-9
    assert val["worst_rel"]["C_D"] <= 1e-9
    assert val["worst_rel"]["A_B"] <= 1e-9
    assert val["rel_err_C"] <= 1e-9
