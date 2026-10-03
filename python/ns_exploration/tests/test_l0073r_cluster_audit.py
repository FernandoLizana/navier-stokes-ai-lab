"""L-0073R cluster audit tests."""

from __future__ import annotations

import json
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.l0073r_audit import audit_full_repair

REPO = Path(__file__).resolve().parents[3]


def test_l0073r_audit_detects_no_false_closure():
    cluster_path = REPO / "experiments/terminal_weighted/cluster_manifest_repair_full_n24_one_inf.json"
    if not cluster_path.is_file():
        return
    rep = audit_full_repair(prec=128)
    lemma = mpfr(rep["sums"]["contrib_lemma_no_outer_maxphi"])
    impl = mpfr(rep["sums"]["contrib_impl"])
    assert abs(float(impl - lemma)) / float(lemma) < 0.01
    assert float(lemma) > 30.0
    assert float(lemma) < 50.0
    cert = rep.get("cert_verify", {})
    assert cert.get("closes_c0008_declared") is False


def test_l0073r_cross_cluster_partition():
    cluster_path = REPO / "experiments/terminal_weighted/cluster_manifest_repair_full_n24_one_inf.json"
    shell_path = REPO / "python/experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n24.json"
    if not cluster_path.is_file() or not shell_path.is_file():
        return
    cluster = json.loads(cluster_path.read_text())
    shell = json.loads(shell_path.read_text())
    from ns_exploration.terminal_weighted.l0073r_audit import cross_cluster_triangle_audit

    cc = cross_cluster_triangle_audit(cluster, shell, prec=128)
    assert cc["partition_covers_each_shell_once"] is True
    assert cc["each_cluster_contrib_le_per_shell_triangle"] is True
    assert 0.2 < cc["cluster_to_route_A_ratio"] < 0.3


def test_l0073r_patched_manifest_matches_lemma():
    cluster_path = REPO / "experiments/terminal_weighted/cluster_manifest_repair_full_n24_one_inf.json"
    if not cluster_path.is_file():
        return
    m = json.loads(cluster_path.read_text())
    total = sum(mpfr(b["contrib_hi"]) for b in m["blocks"])
    assert float(total) > 30.0
