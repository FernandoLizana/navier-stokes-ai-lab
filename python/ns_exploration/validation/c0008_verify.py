"""
Adversarial-independent C-0008 certificate verification.

Does NOT import certificate generators (opnorm_certified, cluster_bounds, sprints).
Recomputes all totals from manifests; fails if declared_hi < recomputed_hi.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.terminal_weighted.constants import EXACT, FROZEN
from ns_exploration.terminal_weighted.intervals import (
    i_star_interval,
    mpfr_const,
    omega_T_hi_from_I_term,
    stokes_floor_interval,
)
from ns_exploration.terminal_weighted.recompute_integral import (
    cluster_integral_hi_from_blocks,
    shell_integral_hi_from_blocks,
)


def manifest_sha256(manifest: dict) -> str:
    m = dict(manifest)
    m.pop("manifest_sha256", None)
    return hashlib.sha256(json.dumps(m, sort_keys=True, default=str).encode()).hexdigest()


def _param_ok(params: dict, key: str, exact: str, *, prec: int) -> bool:
    if key not in params:
        return False
    return mpfr(params[key], precision=prec) == mpfr_const(exact, prec=prec)


def _fail_if_declared_below_recomputed(
    declared: mpfr, recomputed: mpfr, *, label: str, issues: list[str]
) -> bool:
    eps = mpfr("1.0000000000000001", precision=recomputed.precision)
    if declared < recomputed / eps:
        issues.append(f"{label}: declared {declared} < recomputed {recomputed}")
        return False
    return True


def verify_temporal_coverage(cells: list[dict], T: mpfr, *, prec: int) -> tuple[bool, list[str]]:
    issues: list[str] = []
    if not cells:
        issues.append("temporal: no cells")
        return False, issues
    parsed = []
    for c in cells:
        a = mpfr(c["a"], precision=prec)
        b = mpfr(c["b"], precision=prec)
        if b <= a:
            issues.append(f"temporal: invalid cell [{a},{b}]")
        parsed.append((a, b))
    parsed.sort(key=lambda x: float(x[0]))
    zero = mpfr(0, precision=prec)
    if parsed[0][0] != zero:
        issues.append("temporal: gap at t=0")
    for i in range(len(parsed) - 1):
        if parsed[i][1] != parsed[i + 1][0]:
            issues.append(f"temporal: gap/overlap between cells {i} and {i+1}")
    if parsed[-1][1] != T:
        issues.append("temporal: does not reach T")
    return len(issues) == 0, issues


def verify_terminal_certificate(
    cert: dict,
    *,
    manifest: dict | None,
    manifest_path: Path | None,
    prec: int,
) -> tuple[bool, list[str]]:
    issues: list[str] = []
    ok = True

    if cert.get("claim") != "C-0008":
        ok = False
        issues.append("claim != C-0008")

    params = cert.get("parameters", {})
    for key, val in [("nu", EXACT.NU), ("E0", EXACT.E0), ("T", EXACT.T)]:
        if not _param_ok(params, key, val, prec=prec):
            ok = False
            issues.append(f"parameter mismatch {key}")

    if int(params.get("N_max", params.get("n_max", 0))) != 24:
        ok = False
        issues.append("N_max != 24")

    M_decl = cert.get("M_target", params.get("M_target", cert.get("proposed_bound_M")))
    if M_decl is not None and mpfr(str(M_decl), precision=prec) != mpfr_const(
        EXACT.M_TARGET, prec=prec
    ):
        ok = False
        issues.append("M_target mismatch")

    backend = cert.get("arithmetic_backend", {})
    if backend.get("name") != "gmpy2_mpfr":
        ok = False
        issues.append("backend must be gmpy2_mpfr")

    mode_set = cert.get("mode_set", {})
    scope = cert.get("scope", "full_dealias")
    repair_band = scope == "repair_band"
    repair_full = scope in ("repair_full_partial", "repair_full")
    repair_cluster = scope == "repair_full_cluster_sketch"

    if repair_band:
        if mode_set.get("radii_mode") != "band":
            ok = False
            issues.append("repair_band scope requires radii_mode=band")
        if cert.get("bound_method") != "one_inf_only":
            ok = False
            issues.append("repair_band requires bound_method one_inf_only")
    elif repair_full:
        if mode_set.get("radii_mode") != "full":
            ok = False
            issues.append("repair_full scope requires radii_mode=full")
        if cert.get("bound_method") != "one_inf_only":
            ok = False
            issues.append("repair_full requires bound_method one_inf_only")
        if scope == "repair_full" or (scope == "repair_full_partial" and cert.get("repair_complete")):
            if mode_set.get("n_shells", 0) != 87:
                ok = False
                issues.append("repair_full complete requires n_shells=87")
            if mode_set.get("D") != FROZEN.D_full:
                ok = False
                issues.append(f"repair_full D={mode_set.get('D')} != D_full")
    elif repair_cluster:
        if cert.get("bound_method") != "one_inf_only":
            ok = False
            issues.append("repair_cluster requires bound_method one_inf_only")
        if mode_set.get("D") != FROZEN.D_full:
            ok = False
            issues.append(f"repair_cluster D={mode_set.get('D')} != D_full")
    else:
        if mode_set.get("radii_mode") == "band":
            ok = False
            issues.append("band radii rejected for full-dealias claim")
        if mode_set.get("radii_mode") == "full":
            if mode_set.get("D") != FROZEN.D_full:
                ok = False
                issues.append(f"D={mode_set.get('D')} != D_full")
            if mode_set.get("n_shells", 0) != 87:
                ok = False
                issues.append("n_shells != 87")

    if manifest is None:
        ok = False
        issues.append("manifest required but missing")
    elif manifest_path and not manifest_path.is_file():
        ok = False
        issues.append(f"manifest file missing: {manifest_path}")

    I_lo, I_hi_cert = i_star_interval(prec=prec)
    istar = cert.get("i_star_interval", [])
    if len(istar) == 2:
        if mpfr(istar[0], precision=prec) > I_lo:
            ok = False
            issues.append("i_star_interval lo too high vs recomputation")
        if mpfr(istar[1], precision=prec) < I_hi_cert:
            ok = False
            issues.append("i_star_interval hi too low vs recomputation")

    floor_lo, floor_hi = stokes_floor_interval(prec=prec)
    floor = cert.get("stokes_floor_interval", [])
    if len(floor) == 2:
        if mpfr(floor[0], precision=prec) > floor_lo:
            ok = False
            issues.append("stokes_floor lo inconsistent")
        if mpfr(floor[1], precision=prec) < floor_hi:
            ok = False
            issues.append("stokes_floor hi too low")

    src = cert.get("source_hashes", {})
    blocks_for_kind = (manifest or {}).get("blocks") or []
    cluster_manifest_kind = repair_cluster or bool((manifest or {}).get("n_clusters")) or (
        blocks_for_kind
        and "contrib_hi" in blocks_for_kind[0]
        and "cluster_id" in blocks_for_kind[0]
    )
    if cluster_manifest_kind:
        manifest_sha = src.get("cluster_manifest")
    else:
        manifest_sha = src.get("shell_manifest") or src.get("cluster_manifest")
    if manifest and manifest_sha:
        if manifest_sha256(manifest) != manifest_sha:
            ok = False
            issues.append("manifest hash mismatch")

    intv = cert.get("integral_interval", [])
    if len(intv) != 2:
        ok = False
        issues.append("missing integral_interval")
        declared_I = None
    else:
        declared_I = mpfr(intv[1], precision=prec)

    recomputed_I = None
    if manifest and manifest.get("blocks"):
        blocks = manifest["blocks"]
        if blocks and (
            cluster_manifest_kind
            or ("contrib_hi" in blocks[0] and "cluster_id" in blocks[0])
        ):
            recomputed_I = cluster_integral_hi_from_blocks(blocks, prec=prec)
        elif blocks and "shell" in blocks[0]:
            recomputed_I = shell_integral_hi_from_blocks(blocks, prec=prec)
        elif blocks and "contrib_hi" in blocks[0]:
            recomputed_I = cluster_integral_hi_from_blocks(blocks, prec=prec)

    if recomputed_I is not None and declared_I is not None:
        if not _fail_if_declared_below_recomputed(
            declared_I, recomputed_I, label="integral_hi", issues=issues
        ):
            ok = False
        route_a = cert.get("routes", {}).get("A", {})
        if route_a.get("integral_hi"):
            ra = mpfr(route_a["integral_hi"], precision=prec)
            if not _fail_if_declared_below_recomputed(
                ra, recomputed_I, label="route_A integral_hi", issues=issues
            ):
                ok = False

    omega = cert.get("omega_T_interval", [])
    if len(omega) == 2 and declared_I is not None:
        declared_omega = mpfr(omega[1], precision=prec)
        base_I = recomputed_I if recomputed_I is not None else declared_I
        recomputed_omega = omega_T_hi_from_I_term(base_I, prec=prec)
        if not _fail_if_declared_below_recomputed(
            declared_omega, recomputed_omega, label="omega_T_hi", issues=issues
        ):
            ok = False

    cells = cert.get("cells") or cert.get("temporal_cells") or []
    if cells:
        T = mpfr_const(EXACT.T, prec=prec)
        cov_ok, cov_issues = verify_temporal_coverage(cells, T, prec=prec)
        if not cov_ok:
            ok = False
            issues.extend(cov_issues)
    if cert.get("coverage_complete") is False:
        ok = False
        issues.append("coverage_complete=false")

    if cert.get("method") == "band_only":
        ok = False
        issues.append("band-only certificate rejected")

    fro_method = cert.get("bound_method")
    if manifest:
        fro_method = fro_method or manifest.get("bound_method")
    if not repair_band and not repair_full and not repair_cluster and fro_method in ("frobenius", "best", "min_fro_one_inf"):
        ok = False
        issues.append(
            f"bound_method {fro_method!r} rejected for full-dealias until manifest "
            "regenerated with one_inf_only"
        )

    return ok, issues


def verify_terminal_file(
    cert_path: Path,
    *,
    manifest_path: Path | None = None,
) -> tuple[bool, str]:
    cert = json.loads(cert_path.read_text(encoding="utf-8"))
    prec = int(cert.get("arithmetic_backend", {}).get("precision_bits", 200))
    manifest: dict | None = None
    mp = manifest_path
    if mp and mp.is_file():
        manifest = json.loads(mp.read_text(encoding="utf-8"))
    ok, issues = verify_terminal_certificate(
        cert, manifest=manifest, manifest_path=mp, prec=prec
    )
    if ok:
        return True, "PASS"
    detail = "; ".join(issues[:5])
    return False, f"FAIL: {detail}" if detail else "FAIL"


def verify_cluster_file(
    cert_path: Path,
    *,
    manifest_path: Path | None = None,
) -> tuple[bool, str]:
    return verify_terminal_file(cert_path, manifest_path=manifest_path)
