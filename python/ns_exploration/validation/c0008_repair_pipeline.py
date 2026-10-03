"""
C-0008 repair pipeline: manifest enrichment, Route A / A', certificate finalize, verify.
"""

from __future__ import annotations

import json
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.adaptive_time_certificate import (
    adaptive_lipschitz_certificate,
)
from ns_exploration.terminal_weighted.certify_shell_integral import (
    shell_integral_certificate,
)
from ns_exploration.terminal_weighted.cluster_bounds import (
    build_cluster_manifest,
    cluster_integral_certificate,
)
from ns_exploration.terminal_weighted.cluster_residual_prototype import (
    run_band6_cluster_residual,
)
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.l0048_audit import l0048_value_classification
from ns_exploration.terminal_weighted.n_coverage import (
    band_majorant_ladder,
    n_domination_analysis,
)
from ns_exploration.terminal_weighted.repair_full_manifest import (
    assemble_repair_manifest,
    save_repair_manifest,
)
from ns_exploration.validation.c0008_repair_band_cert import _build_repair_certificate
from ns_exploration.validation.c0008_verify import manifest_sha256, verify_terminal_certificate


def enrich_manifest_metadata(manifest: dict, *, prec: int = 128) -> dict:
    """Attach n_coverage, band ladder, optional band-6 cluster audit."""
    n = int(manifest.get("n", 24))
    out = dict(manifest)
    out["n_coverage"] = n_domination_analysis(n, run_band_ladder=n >= 24)
    if manifest.get("bound_method") == "one_inf_only":
        try:
            out["band_majorant_ladder_16"] = band_majorant_ladder((1, 2, 3, 4, 5, 6), prec=prec)
        except Exception as exc:
            out["band_majorant_ladder_16"] = {"error": str(exc)}
    if len(manifest.get("blocks", [])) >= 6 and manifest.get("n_shells", 0) <= 6:
        try:
            cr = run_band6_cluster_residual(n=n, prec=prec)
            out["cluster_residual_L0073R"] = cr["residual_combined"]
            out["B_C_combined_hi"] = cr["B_C_combined_hi"]
        except Exception as exc:
            out["cluster_residual_L0073R"] = {"error": str(exc)}
    payload = dict(out)
    payload.pop("manifest_sha256", None)
    import hashlib

    out["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, default=str).encode()
    ).hexdigest()
    return out


def build_cluster_manifest_repair(
    per_shell_manifest: dict,
    *,
    prec: int = 128,
    n_clusters: int = 12,
    workers: int = 1,
) -> dict:
    n = int(per_shell_manifest["n"])
    radii = tuple(all_dealias_radii(n))
    if per_shell_manifest.get("n_shells", 0) < len(radii):
        radii = tuple(int(b["shell"]) for b in per_shell_manifest["blocks"])
    return build_cluster_manifest(
        n,
        radii=radii,
        n_clusters=min(n_clusters, max(1, len(radii) // 2)),
        per_shell_manifest=per_shell_manifest,
        bound_method="one_inf_only",
        prec=prec,
        workers=workers,
    )


def finalize_repair_certificate(
    per_shell_manifest: dict,
    *,
    prec: int = 128,
    cluster_manifest: dict | None = None,
    workers: int = 1,
) -> dict:
    """Build repair certificate with Route A, optional A', B sketch."""
    n = int(per_shell_manifest.get("n", 24))
    complete = bool(per_shell_manifest.get("repair_complete"))
    n_total = len(all_dealias_radii(n))
    n_done = len(per_shell_manifest.get("blocks", []))
    is_full = n_done >= n_total or per_shell_manifest.get("n_shells", 0) >= n_total

    if is_full and complete:
        scope = "repair_full"
    elif per_shell_manifest.get("n_shells", 0) <= 6 and n_done <= 6:
        scope = "repair_band"
    else:
        scope = "repair_full_partial"

    route_a = shell_integral_certificate(n, prec=prec, manifest=per_shell_manifest)
    routes: dict = {"A": route_a}

    if cluster_manifest is None and n_done >= 3:
        try:
            cluster_manifest = build_cluster_manifest_repair(
                per_shell_manifest, prec=prec, workers=workers
            )
        except Exception:
            cluster_manifest = None

    if cluster_manifest is not None:
        route_ap = cluster_integral_certificate(cluster_manifest, prec=prec)
        routes["A_prime"] = route_ap

    try:
        if not (is_full and complete):
            route_b = adaptive_lipschitz_certificate(
                n, prec=prec, manifest=per_shell_manifest, initial_cells=8
            )
            routes["B"] = route_b
    except Exception:
        pass

    best_name = min(routes, key=lambda k: mpfr(routes[k]["integral_hi"]))
    best = routes[best_name]

    cert = _build_repair_certificate(
        per_shell_manifest,
        prec=prec,
        scope=scope,
        radii=tuple(range(1, 7)) if scope == "repair_band" else None,
    )
    cert["scope"] = scope
    cert["repair_complete"] = complete and is_full
    cert["routes"] = {
        k: {
            "method": v.get("method"),
            "integral_hi": v.get("integral_hi"),
            "omega_T_hi": v.get("omega_T_hi"),
            "closes_strict": v.get("closes_strict"),
        }
        for k, v in routes.items()
    }
    cert["best_route"] = best_name
    cert["best_route_sketch"] = {
        "integral_hi": best["integral_hi"],
        "omega_T_hi": best.get("omega_T_hi"),
        "note": "Sketch only — declared integral_interval is Route A shell sum (verifier-recomputable).",
    }
    # Declared intervals stay Route A from _build_repair_certificate (adversarial recompute).
    cert["n_coverage"] = n_domination_analysis(n, run_band_ladder=True)
    cert["l0048_audit"] = l0048_value_classification()
    if scope == "repair_band":
        cert["coverage_complete"] = True
    else:
        cert["coverage_complete"] = complete and is_full
    if cluster_manifest is not None:
        cert["source_hashes"]["cluster_manifest"] = cluster_manifest.get("manifest_sha256")
    cert["closes_c0008"] = bool(
        complete
        and is_full
        and mpfr(cert["integral_interval"][1]) < mpfr(cert["i_star_interval"][0])
        and mpfr(cert["omega_T_interval"][1]) < mpfr(FROZEN.c0008_target_M)
    )
    cert["limitations"] = [
        "Repair certificate: one_inf_only per-shell bounds; L-0073R on Route A'.",
        "Historical L-0072/L-0073 manifests are NOT superseded until adversarial PASS.",
        "N<=24 inclusion: N7_sketch (L-0076 empirical ladder).",
    ]
    if not complete:
        cert["limitations"].insert(
            0, f"Partial: {n_done}/{n_total} shells — NOT full-dealias closure."
        )
    return cert


def run_finalize_pipeline(
    manifest_path: Path,
    *,
    prec: int = 128,
    workers: int = 1,
    out_cert: Path | None = None,
    out_cluster: Path | None = None,
) -> dict:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest = enrich_manifest_metadata(manifest, prec=prec)
    save_repair_manifest(manifest, manifest_path)

    cluster_manifest = None
    if len(manifest.get("blocks", [])) >= 3:
        cluster_manifest = build_cluster_manifest_repair(
            manifest, prec=prec, workers=workers
        )
        if out_cluster:
            out_cluster.parent.mkdir(parents=True, exist_ok=True)
            out_cluster.write_text(json.dumps(cluster_manifest, indent=2), encoding="utf-8")

    cert = finalize_repair_certificate(
        manifest, prec=prec, cluster_manifest=cluster_manifest, workers=workers
    )
    n = int(manifest.get("n", 24))
    cert_path = out_cert or Path(
        f"certificates/CERT-C0008-repair-full-n{n}-one-inf.json"
    )
    cert_path.parent.mkdir(parents=True, exist_ok=True)
    cert_path.write_text(json.dumps(cert, indent=2), encoding="utf-8")

    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=manifest_path, prec=prec
    )
    return {
        "manifest": str(manifest_path),
        "cluster_manifest": str(out_cluster) if out_cluster else None,
        "cert": str(cert_path),
        "verify_ok": ok,
        "issues": issues,
        "n_blocks": len(manifest.get("blocks", [])),
        "repair_complete": manifest.get("repair_complete"),
        "best_route": cert.get("best_route"),
        "integral_hi": cert["integral_interval"][1],
        "closes_c0008": cert.get("closes_c0008"),
    }
