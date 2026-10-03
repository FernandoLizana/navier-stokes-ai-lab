"""Phase D certification tests (spec §12)."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.certify_shell_integral import shell_integral_certificate
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval
from ns_exploration.terminal_weighted.l0048_audit import l0048_value_classification
from ns_exploration.terminal_weighted.n_coverage import n_domination_analysis
from ns_exploration.terminal_weighted.opnorm_certified import certified_L_op_frobenius_hi
from ns_exploration.terminal_weighted.shell_decomposition import (
    audit_shell_reconstruction_band,
    build_shell_manifest,
)


@pytest.fixture(scope="module")
def band_radii():
    return tuple(range(1, 7))


def test_l0048_not_rigorous_upper():
    audit = l0048_value_classification(validate_dense=False)
    assert audit["classification"] == "a_approximate_float64_operator_norm"
    assert "c_rigorous_upper" in audit["not_certified_as"]


def test_single_pass_matches_per_shell(band_radii):
    from ns_exploration.terminal_weighted.opnorm_certified import (
        certified_L_op_frobenius_hi,
        certified_shell_blocks,
    )

    b1 = certified_L_op_frobenius_hi(24, band_radii, output_shell=1, prec=128)
    blk = certified_shell_blocks(24, band_radii, prec=128, bound_method="frobenius")
    b1sp = next(b for b in blk["blocks"] if b["shell"] == 1)
    assert b1.C_term_hi == b1sp["C_term_hi"]


def test_recompute_integral_from_blocks(band_radii):
    from ns_exploration.terminal_weighted.recompute_integral import (
        shell_integral_hi_from_blocks,
    )

    manifest = build_shell_manifest(24, prec=128, radii=band_radii)
    I_re = shell_integral_hi_from_blocks(manifest["blocks"], prec=128)
    cert = shell_integral_certificate(24, prec=128, manifest=manifest)
    assert mpfr(cert["integral_hi"]) <= I_re * mpfr("1.0000000000000001")


def test_i_star_recomputed():
    I_lo, I_hi = i_star_interval(prec=200)
    assert I_lo < I_hi
    assert float(I_hi) < FROZEN.integral_threshold * 1.001
    assert float(I_lo) > FROZEN.integral_threshold * 0.999


def test_n_coverage_wavevectors_subset():
    cov = n_domination_analysis(24)
    assert cov["all_subset"]
    assert cov["evidence_level"] == "N7_sketch"


def test_shell_reconstruction_band(band_radii):
    rec = audit_shell_reconstruction_band(24, band_radii)
    assert rec["direct_vs_tensor_ok"]
    assert rec["frobenius_full"] > 0


def test_certified_frobenius_single_shell(band_radii):
    bnd = certified_L_op_frobenius_hi(24, band_radii, output_shell=1, prec=128)
    assert float(mpfr(bnd.C_term_hi)) > 0
    assert bnd.evidence_level == "N5"
    assert bnd.method == "frobenius_mpfr_up"


def test_manifest_hash_mutates(band_radii):
    m1 = build_shell_manifest(24, prec=128, radii=band_radii)
    m2 = dict(m1)
    m2["blocks"] = list(m2["blocks"])
    m2["blocks"][0] = dict(m2["blocks"][0])
    m2["blocks"][0]["C_term_hi"] = "999"
    import hashlib

    d2 = hashlib.sha256(json.dumps(m2, sort_keys=True, default=str).encode()).hexdigest()
    assert d2 != m1["manifest_sha256"]


def test_shell_integral_band_pilot(band_radii):
    manifest = build_shell_manifest(24, prec=128, radii=band_radii)
    cert = shell_integral_certificate(24, prec=128, manifest=manifest)
    assert cert["method"] == "shell_integral_route_A"
    assert float(mpfr(cert["integral_hi"])) > 0
    assert cert["full_dealias"] is False
    assert cert["closes_c0008"] is False


def test_verifier_rejects_open_integral(tmp_path, band_radii):
    manifest = build_shell_manifest(24, prec=128, radii=band_radii)
    cert = shell_integral_certificate(24, prec=128, manifest=manifest)
    I_lo, I_hi = i_star_interval(prec=128)
    payload = {
        "claim": "C-0008",
        "status": "candidate_certificate",
        "parameters": {"nu": FROZEN.nu, "E0": FROZEN.E0, "T": FROZEN.T},
        "i_star_interval": [str(I_lo), str(I_hi)],
        "stokes_floor_interval": [FROZEN.stokes_floor, FROZEN.stokes_floor + 1e-9],
        "integral_interval": ["0", cert["integral_hi"]],
        "omega_T_interval": ["0", cert["omega_T_hi"]],
        "method": "shell_integral_route_A",
    }
    p = tmp_path / "cert.json"
    p.write_text(json.dumps(payload), encoding="utf-8")

    import sys

    root = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(root))
    from verify_c0008_terminal_certificate import verify_certificate

    ok, msgs = verify_certificate(p)
    assert not ok
    assert any("FAIL" in m for m in msgs)


def test_verifier_rejects_band_only_method(tmp_path):
    I_lo, I_hi = i_star_interval()
    payload = {
        "claim": "C-0008",
        "status": "candidate_certificate",
        "parameters": {"nu": FROZEN.nu, "E0": FROZEN.E0, "T": FROZEN.T},
        "i_star_interval": [str(I_lo), str(I_hi)],
        "stokes_floor_interval": [FROZEN.stokes_floor, FROZEN.stokes_floor],
        "integral_interval": ["0", "0.5"],
        "omega_T_interval": ["0", "41.0"],
        "method": "band_only",
    }
    p = tmp_path / "cert.json"
    p.write_text(json.dumps(payload), encoding="utf-8")
    import sys

    root = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(root))
    from verify_c0008_terminal_certificate import verify_certificate

    ok, msgs = verify_certificate(p)
    assert not ok
    assert any("band-only" in m for m in msgs)


def test_one_inf_tighter_than_frobenius_band(band_radii):
    manifest = build_shell_manifest(24, prec=128, radii=band_radii, workers=1, bound_method="best")
    for blk in manifest["blocks"]:
        assert float(mpfr(blk["C_term_one_inf_hi"])) <= float(mpfr(blk["C_term_frobenius_hi"]))
    assert float(mpfr(manifest["full_best_hi"]["C_term_hi"])) <= float(
        mpfr(manifest["full_frobenius_hi"]["C_term_hi"])
    )


def test_one_inf_single_shell_matches_band(band_radii):
    from ns_exploration.terminal_weighted.opnorm_certified import (
        certified_L_op_one_inf_single_shell_hi,
    )

    b = certified_L_op_one_inf_single_shell_hi(24, band_radii, 1, prec=128)
    manifest = build_shell_manifest(24, prec=128, radii=band_radii, bound_method="best")
    ref = next(x for x in manifest["blocks"] if x["shell"] == 1)
    assert b.C_term_hi == ref["C_term_one_inf_hi"]


def test_cluster_frobenius_tighter_than_shell_sum(band_radii):
    from ns_exploration.terminal_weighted.opnorm_certified import compare_shell_vs_cluster_band

    cmp = compare_shell_vs_cluster_band(24, band_radii, prec=128)
    assert float(cmp["cluster_tighter_ratio"]) <= 1.0 + 1e-9
    assert float(cmp["cluster_C_term_hi"]) > 0


def test_route_a_prime_band_improves_integral(band_radii):
    from ns_exploration.terminal_weighted.cluster_bounds import (
        build_cluster_manifest,
        cluster_integral_certificate,
    )

    per_shell = build_shell_manifest(24, prec=128, radii=band_radii, bound_method="frobenius")
    manifest = build_cluster_manifest(
        24,
        radii=band_radii,
        n_clusters=2,
        prec=128,
        workers=1,
        per_shell_manifest=per_shell,
    )
    cert = cluster_integral_certificate(manifest, prec=128)
    from ns_exploration.terminal_weighted.recompute_integral import shell_integral_hi_from_blocks

    I_shell = shell_integral_hi_from_blocks(per_shell["blocks"], prec=128)
    I_cluster = mpfr(cert["integral_hi"])
    assert I_cluster <= I_shell * mpfr("1.0000000000000001")


def test_cluster_route_d_band(band_radii):
    from ns_exploration.terminal_weighted.cluster_bounds import (
        build_cluster_manifest,
        cluster_route_d_certificate,
    )

    per_shell = build_shell_manifest(24, prec=128, radii=band_radii, bound_method="best")
    manifest = build_cluster_manifest(
        24,
        radii=band_radii,
        n_clusters=2,
        prec=128,
        workers=1,
        per_shell_manifest=per_shell,
        bound_method="best",
    )
    cert_d = cluster_route_d_certificate(manifest, per_shell_manifest=per_shell, prec=128)
    assert float(cert_d["integral_hi"]) > 0
    assert float(cert_d["omega_T_hi"]) > float(FROZEN.stokes_floor)


def test_adaptive_partitions_band(band_radii):
    from ns_exploration.terminal_weighted.cluster_bounds import adaptive_cluster_partition

    for strategy in ("equal", "fine_low", "singleton"):
        parts = adaptive_cluster_partition(band_radii, strategy, n_clusters=2)
        covered = sorted(r for p in parts for r in p.shells)
        assert covered == sorted(set(band_radii))


def test_seed_equal_12_from_l0073():
    from ns_exploration.experiments.sprint_c0008_route_a_ladder import seed_equal_12_from_l0073

    row = seed_equal_12_from_l0073(prec=128)
    if row is None:
        pytest.skip("L-0073 manifest missing")
    assert row["label"] == "equal_12"
    assert float(row["route_a_I_hi"]) < float(row["route_d_I_hi"])
    assert row["best_route"] == "A"


def test_verifier_parameter_nu_fail(tmp_path):
    I_lo, I_hi = i_star_interval()
    payload = {
        "claim": "C-0008",
        "status": "candidate_certificate",
        "parameters": {"nu": 0.2, "E0": FROZEN.E0, "T": FROZEN.T, "N_max": 24},
        "i_star_interval": [str(I_lo), str(I_hi)],
        "stokes_floor_interval": [FROZEN.stokes_floor, FROZEN.stokes_floor],
        "integral_interval": ["0", "0.5"],
        "omega_T_interval": ["0", "41.0"],
        "method": "shell_integral_route_A",
        "mode_set": {"radii_mode": "full", "D": FROZEN.D_full, "n_shells": 87},
        "arithmetic_backend": {"name": "gmpy2_mpfr"},
    }
    p = tmp_path / "cert.json"
    p.write_text(json.dumps(payload), encoding="utf-8")
    import sys
    root = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(root))
    from verify_c0008_terminal_certificate import verify_certificate
    ok, msgs = verify_certificate(p)
    assert not ok
    assert any("parameter mismatch" in m for m in msgs)
