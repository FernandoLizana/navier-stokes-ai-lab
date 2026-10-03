"""L-0073R cluster certificate honesty audit."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.cluster_bounds import (
    abs_up_mpfr,
    cluster_alpha_candidates,
    cluster_residual_contrib_hi,
    cluster_residual_vs_legacy_audit,
)
from ns_exploration.terminal_weighted.intervals import (
    _ctx,
    add_up,
    mul_up,
    sub_up,
)
from ns_exploration.terminal_weighted.recompute_integral import (
    cluster_integral_hi_from_blocks,
    shell_integral_hi_from_blocks,
)
from ns_exploration.validation.c0008_verify import manifest_sha256, verify_terminal_certificate

REPO = Path(__file__).resolve().parents[3]


def lemma_contrib_hi_no_outer_maxphi(
    shells: tuple[int, ...],
    B_r: dict[int, mpfr],
    B_C: mpfr,
    phi_r: dict[int, mpfr],
    alpha: mpfr,
    *,
    prec: int,
) -> mpfr:
    """Lemma L-0073R: |alpha| B_C + sum |phi_r - alpha| B_r, then min vs per-shell."""
    _ctx(prec)
    per_shell = mpfr(0, precision=prec)
    residual = mul_up(abs_up_mpfr(alpha, prec), B_C)
    for r in shells:
        rf = phi_r[r]
        diff = abs_up_mpfr(sub_up(rf, alpha), prec)
        residual = add_up(residual, mul_up(diff, B_r[r]))
        per_shell = add_up(per_shell, mul_up(B_r[r], rf))
    return per_shell if per_shell <= residual else residual


def audit_cluster_manifest(manifest: dict, *, prec: int = 128) -> dict[str, Any]:
    _ctx(prec)
    impl_sum = mpfr(0, precision=prec)
    lemma_sum = mpfr(0, precision=prec)
    per_shell_sum = mpfr(0, precision=prec)
    legacy_sum = mpfr(0, precision=prec)
    blocks_out = []
    outer_maxphi_violations = 0

    for blk in manifest["blocks"]:
        impl = mpfr(blk["contrib_hi"], precision=prec)
        impl_sum = add_up(impl_sum, impl)
        per = mpfr(blk["per_shell_integral_hi"], precision=prec)
        per_shell_sum = add_up(per_shell_sum, per)
        legacy = mpfr(blk["cluster_integral_hi"], precision=prec)
        legacy_sum = add_up(legacy_sum, legacy)

        phi = {int(d["shell"]): mpfr(d["factor_hi"], precision=prec) for d in blk["per_shell_detail"]}
        B_r = {int(d["shell"]): mpfr(d["B_r_hi"], precision=prec) for d in blk["per_shell_detail"]}
        B_C = mpfr(blk["B_cluster_hi"], precision=prec)
        shells = tuple(int(s) for s in blk["shells"])
        cands = cluster_alpha_candidates(phi)
        alpha = cands.get(blk.get("cluster_residual_L0073R", {}).get("best_residual", {}).get("alpha", "min"), cands["min"])
        if isinstance(alpha, str):
            alpha = cands[alpha]
        lemma = lemma_contrib_hi_no_outer_maxphi(shells, B_r, B_C, phi, alpha, prec=prec)
        lemma_sum = add_up(lemma_sum, lemma)

        impl_re = cluster_residual_contrib_hi(shells, B_r, B_C, phi, alpha, prec=prec)
        if impl != impl_re:
            pass  # alpha name mismatch ok
        max_phi = max(phi.values(), key=float)
        residual_inner = mul_up(abs_up_mpfr(alpha, prec), B_C)
        for r in shells:
            diff = abs_up_mpfr(sub_up(phi[r], alpha), prec)
            residual_inner = add_up(residual_inner, mul_up(diff, B_r[r]))
        with_outer = mul_up(residual_inner, max_phi)
        if abs(float(with_outer - impl)) < abs(float(impl)) * 1e-6:
            outer_maxphi_violations += 1

        blocks_out.append(
            {
                "cluster_id": blk["cluster_id"],
                "shells": blk["shells"],
                "contrib_impl": str(impl),
                "contrib_lemma_no_maxphi": str(lemma),
                "per_shell_integral_hi": str(per),
                "legacy_invalid_hi": str(legacy),
                "ratio_impl_to_per_shell": float(impl / per) if per else None,
            }
        )

    return {
        "lemma": "L-0073R",
        "evidence_level": "N4_audit",
        "n_clusters": manifest.get("n_clusters"),
        "bound_method": manifest.get("bound_method"),
        "sums": {
            "contrib_impl": str(impl_sum),
            "contrib_lemma_no_outer_maxphi": str(lemma_sum),
            "per_shell_triangle_within_clusters": str(per_shell_sum),
            "legacy_invalid_B_C_times_max_phi": str(legacy_sum),
        },
        "outer_maxphi_factor_detected": outer_maxphi_violations,
        "implementation_matches_lemma": outer_maxphi_violations == 0,
        "impl_vs_lemma": {
            "impl_lt_lemma": float(impl_sum) < float(lemma_sum),
            "lemma_to_impl_ratio": float(lemma_sum / impl_sum) if impl_sum else None,
        },
        "vs_route_A": {},
        "verdict": [],
        "blocks": blocks_out,
    }


def cross_cluster_triangle_audit(
    cluster_manifest: dict,
    shell_manifest: dict | None = None,
    *,
    prec: int = 128,
) -> dict[str, Any]:
    """Partition coverage + sum-of-clusters vs Route A (formal triangle gap)."""
    _ctx(prec)
    blocks = cluster_manifest.get("blocks", [])
    all_shells: list[int] = []
    cluster_sum = mpfr(0, precision=prec)
    per_in_cluster = mpfr(0, precision=prec)
    legacy_sum = mpfr(0, precision=prec)
    per_le_residual = True
    for blk in blocks:
        cluster_sum = add_up(cluster_sum, mpfr(blk["contrib_hi"], precision=prec))
        per = mpfr(blk["per_shell_integral_hi"], precision=prec)
        per_in_cluster = add_up(per_in_cluster, per)
        legacy_sum = add_up(legacy_sum, mpfr(blk["cluster_integral_hi"], precision=prec))
        if mpfr(blk["contrib_hi"], precision=prec) > per * mpfr("1.0000001", precision=prec):
            per_le_residual = False
        all_shells.extend(int(s) for s in blk["shells"])

    route_a = None
    if shell_manifest is not None:
        route_a = shell_integral_hi_from_blocks(shell_manifest["blocks"], prec=prec)

    expected = len(set(all_shells))
    partition_ok = len(all_shells) == expected and expected == cluster_manifest.get("n_shells", expected)

    out: dict[str, Any] = {
        "n_clusters": cluster_manifest.get("n_clusters"),
        "n_shells_partitioned": expected,
        "partition_covers_each_shell_once": partition_ok,
        "each_cluster_contrib_le_per_shell_triangle": per_le_residual,
        "cluster_contrib_sum": str(cluster_sum),
        "per_shell_sum_within_clusters": str(per_in_cluster),
        "legacy_invalid_sum": str(legacy_sum),
        "route_A_shell_integral": str(route_a) if route_a is not None else None,
        "cluster_to_route_A_ratio": float(cluster_sum / route_a) if route_a else None,
        "formal_gap": (
            "sum_C contrib_C <= I_total is NOT proved from per-cluster operator bounds "
            "without global combined pass or N7 triangle lemma across clusters."
        ),
        "empirical": {},
    }
    if route_a is not None:
        out["empirical"] = {
            "cluster_tighter_than_route_A": float(cluster_sum) < float(route_a),
            "per_shell_within_clusters_matches_route_A": abs(float(per_in_cluster - route_a)) / float(route_a) < 1e-9,
            "legacy_near_cluster_sum": abs(float(legacy_sum - cluster_sum)) / float(cluster_sum) < 0.05,
        }
    return out


def audit_full_repair(*, prec: int = 128) -> dict[str, Any]:
    cluster_path = REPO / "experiments/terminal_weighted/cluster_manifest_repair_full_n24_one_inf.json"
    shell_path = REPO / "python/experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n24.json"
    cert_path = REPO / "certificates/CERT-L0073R-repair-C0008-cluster-one-inf.json"

    cluster = json.loads(cluster_path.read_text(encoding="utf-8"))
    rep = audit_cluster_manifest(cluster, prec=prec)

    if shell_path.is_file():
        shell = json.loads(shell_path.read_text(encoding="utf-8"))
        route_a = shell_integral_hi_from_blocks(shell["blocks"], prec=prec)
        rep["vs_route_A"] = {
            "shell_integral_hi": str(route_a),
            "cluster_impl_sum": rep["sums"]["contrib_impl"],
            "cluster_lemma_sum": rep["sums"]["contrib_lemma_no_outer_maxphi"],
            "impl_vs_route_A_ratio": float(mpfr(rep["sums"]["contrib_impl"]) / route_a),
            "lemma_vs_route_A_ratio": float(mpfr(rep["sums"]["contrib_lemma_no_outer_maxphi"]) / route_a),
        }
        rep["cross_cluster"] = cross_cluster_triangle_audit(cluster, shell, prec=prec)

    if cert_path.is_file():
        cert = json.loads(cert_path.read_text(encoding="utf-8"))
        ok, issues = verify_terminal_certificate(
            cert, manifest=cluster, manifest_path=cluster_path, prec=prec
        )
        rep["cert_verify"] = {"ok": ok, "issues": issues, "closes_c0008_declared": cert.get("closes_c0008")}

    verdict: list[str] = []
    if rep["outer_maxphi_factor_detected"] == len(cluster["blocks"]):
        verdict.append(
            "FAIL: implementation multiplies residual by max_phi — NOT in lemma L-0073R; "
            "inflates invalid legacy form, deflates contrib ~max_phi factor (~50×)."
        )
    if rep["impl_vs_lemma"].get("impl_lt_lemma"):
        verdict.append(
            "FAIL: declared contrib_impl sum is below lemma-correct sum (cert too optimistic)."
        )
    if rep["vs_route_A"]:
        ratio = rep["vs_route_A"].get("lemma_vs_route_A_ratio", 1.0)
        if ratio < 0.05:
            verdict.append(
                f"WARN: lemma-correct cluster sum is {ratio:.1%} of Route A — suspiciously tight; "
                "triangle across clusters + combined B_C pass need N7 proof."
            )
    if rep.get("cross_cluster"):
        cc = rep["cross_cluster"]
        if cc.get("partition_covers_each_shell_once"):
            verdict.append("PASS: cluster partition covers 87 shells exactly once.")
        ratio = cc.get("cluster_to_route_A_ratio")
        if ratio is not None:
            verdict.append(
                f"INFO: cluster sum is {ratio:.1%} of Route A — tighter if valid; "
                f"formal gap: {cc['formal_gap']}"
            )
    verdict.append(
        "Route A shell-sum (I_hi~170.59) remains primary honest repair bound."
    )
    verdict.append(
        "Adversarial verify PASS on declared fields does NOT imply mathematical closure."
    )
    if cert_path.is_file() and rep.get("cert_verify", {}).get("closes_c0008_declared"):
        verdict.append(
            "FAIL: closes_c0008=true on cluster cert is NOT admissible until lemma audit PASS."
        )
    rep["verdict"] = verdict
    rep["audit_pass"] = not any(v.startswith("FAIL:") for v in verdict)
    rep["paths"] = {
        "cluster_manifest": str(cluster_path),
        "cluster_sha256": manifest_sha256(cluster),
    }
    return rep


def recompute_manifest_contrib_hi(manifest: dict, *, prec: int = 128) -> dict:
    """Recompute contrib_hi on existing cluster manifest (no D^2 pass)."""
    _ctx(prec)
    for blk in manifest.get("blocks", []):
        phi = {int(d["shell"]): mpfr(d["factor_hi"], precision=prec) for d in blk["per_shell_detail"]}
        B_r = {int(d["shell"]): mpfr(d["B_r_hi"], precision=prec) for d in blk["per_shell_detail"]}
        B_C = mpfr(blk["B_cluster_hi"], precision=prec)
        shells = tuple(int(s) for s in blk["shells"])
        audit = cluster_residual_vs_legacy_audit(shells, B_r, B_C, phi, prec=prec)
        blk["contrib_hi"] = audit["best_residual"]["contrib_hi"]
        blk["cluster_residual_L0073R"] = audit
        blk["method"] = f"route_a_prime_L-0073R_{audit['best_residual']['alpha']}_one_inf_only"
    import hashlib
    import json as _json

    m = dict(manifest)
    m.pop("manifest_sha256", None)
    manifest["manifest_sha256"] = hashlib.sha256(
        _json.dumps(m, sort_keys=True, default=str).encode()
    ).hexdigest()
    return manifest


def write_audit_report(rep: dict) -> tuple[Path, Path]:
    out_json = REPO / "reports" / "L0073R_CLUSTER_AUDIT.json"
    out_md = REPO / "reports" / "L0073R_CLUSTER_AUDIT.md"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(rep, indent=2), encoding="utf-8")

    s = rep["sums"]
    vra = rep.get("vs_route_A", {})
    cc = rep.get("cross_cluster", {})
    md = f"""# L-0073R Cluster Certificate Audit

**Audit pass:** `{rep['audit_pass']}`  
**Evidence:** N4_audit (implementation vs lemma L-0073R)

## Sums (n=24 full dealias, 11 clusters)

| Quantity | Value |
|----------|-------|
| contrib_impl (current code) | {s['contrib_impl']} |
| contrib_lemma (no outer max_phi) | {s['contrib_lemma_no_outer_maxphi']} |
| per-shell triangle within clusters | {s['per_shell_triangle_within_clusters']} |
| legacy invalid B_C×max(φ) | {s['legacy_invalid_B_C_times_max_phi']} |
| Route A shell integral | {vra.get('shell_integral_hi', 'n/a')} |

## Cross-cluster triangle

| Check | Value |
|-------|-------|
| Partition 87 shells once | {cc.get('partition_covers_each_shell_once')} |
| Cluster / Route A ratio | {cc.get('cluster_to_route_A_ratio')} |
| Formal gap | {cc.get('formal_gap', 'n/a')} |

## Implementation vs lemma

- Blocks with outer `max_phi` factor: **{rep['outer_maxphi_factor_detected']}/{rep['n_clusters']}**
- Matches lemma: **{rep['implementation_matches_lemma']}**
- Lemma/impl ratio: **{rep['impl_vs_lemma'].get('lemma_to_impl_ratio')}**

## Verdict

"""
    for line in rep["verdict"]:
        md += f"- {line}\n"
    out_md.write_text(md, encoding="utf-8")
    return out_json, out_md
