"""
Cluster partition and Route A' integral bounds for terminal-weighted operator.

Rigorous per-cluster contribution:
  I_C <= min( sum_{r in C} B_r f(r),  B_C * max_{r in C} f(r) )

where B_r from per-shell manifest, B_C from combined-output Frobenius pass,
f(r) = (1 - exp(-2 nu r T)) / (2 nu r).
"""

from __future__ import annotations

import hashlib
import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import (
    ArithmeticBackend,
    add_up,
    div_up,
    exp_down,
    i_star_interval,
    mul_up,
    sqrt_up,
    sub_up,
    _ctx,
)
from ns_exploration.terminal_weighted.opnorm_certified import (
    _accumulate_cluster_bounds_single_pass,
    _accumulate_frobenius_clusters_single_pass,
    _bound_from_fro_sq,
    _bound_from_one_inf_rowlist,
)


@dataclass
class ClusterSpec:
    cluster_id: int
    shells: tuple[int, ...]

    def as_dict(self) -> dict:
        return {"cluster_id": self.cluster_id, "shells": list(self.shells)}


def equal_size_cluster_partition(
    radii: tuple[int, ...], n_clusters: int = 12
) -> list[ClusterSpec]:
    """Split sorted dealias shells into n_clusters contiguous groups by index."""
    shells = sorted(set(radii))
    n_clusters = max(1, min(n_clusters, len(shells)))
    chunk = max(1, (len(shells) + n_clusters - 1) // n_clusters)
    out: list[ClusterSpec] = []
    cid = 0
    for i in range(0, len(shells), chunk):
        group = tuple(shells[i : i + chunk])
        out.append(ClusterSpec(cluster_id=cid, shells=group))
        cid += 1
    return out


def singleton_cluster_partition(radii: tuple[int, ...]) -> list[ClusterSpec]:
    return [
        ClusterSpec(cluster_id=i, shells=(r,))
        for i, r in enumerate(sorted(set(radii)))
    ]


def fine_low_cluster_partition(
    radii: tuple[int, ...], *, low_cutoff: int = 12
) -> list[ClusterSpec]:
    """Each shell r<=low_cutoff alone; higher shells in equal groups of ~8."""
    shells = sorted(set(radii))
    low = [r for r in shells if r <= low_cutoff]
    high = [r for r in shells if r > low_cutoff]
    out: list[ClusterSpec] = []
    cid = 0
    for r in low:
        out.append(ClusterSpec(cluster_id=cid, shells=(r,)))
        cid += 1
    if high:
        chunk = max(1, (len(high) + 7) // 8)
        for i in range(0, len(high), chunk):
            out.append(ClusterSpec(cluster_id=cid, shells=tuple(high[i : i + chunk])))
            cid += 1
    return out


def adaptive_cluster_partition(
    radii: tuple[int, ...],
    strategy: str = "equal",
    *,
    n_clusters: int = 12,
    low_cutoff: int = 12,
) -> list[ClusterSpec]:
    if strategy == "singleton":
        return singleton_cluster_partition(radii)
    if strategy == "fine_low":
        return fine_low_cluster_partition(radii, low_cutoff=low_cutoff)
    if strategy == "equal":
        return equal_size_cluster_partition(radii, n_clusters=n_clusters)
    raise ValueError(f"unknown cluster strategy: {strategy!r}")


def shell_to_cluster_map(clusters: list[ClusterSpec]) -> dict[int, int]:
    m: dict[int, int] = {}
    for c in clusters:
        for r in c.shells:
            m[r] = c.cluster_id
    return m


def _shell_factor_hi(r: mpfr, *, nu: mpfr, T: mpfr, prec: int) -> mpfr:
    two = mpfr(2, precision=prec)
    expo = mul_up(mul_up(mul_up(mpfr(-2, precision=prec), nu), r), T)
    e_down = exp_down(expo)
    with _ctx(prec, up=True):
        factor_num_hi = mpfr(1) - e_down
    denom = mul_up(mul_up(two, nu), r)
    return div_up(factor_num_hi, denom)


def _shell_correction_hi(r: mpfr, *, nu: mpfr, T: mpfr, prec: int) -> mpfr:
    """Route D factor int_phi(r) = T - (1-e^{-2 nu r T})/(2 nu r), upward rounded."""
    from ns_exploration.terminal_weighted.intervals import div_down, exp_up, sub_down, sub_up

    two = mpfr(2, precision=prec)
    expo_neg = mul_up(mul_up(mul_up(two, nu), r), T)
    e_up = exp_up(mul_up(mpfr(-1, precision=prec), expo_neg))
    one = mpfr(1, precision=prec)
    correction_lo = div_down(sub_down(one, e_up), mul_up(mul_up(two, nu), r))
    return sub_up(T, correction_lo)


def cluster_residual_contrib_hi(
    shells: tuple[int, ...],
    B_r: dict[int, mpfr],
    B_C: mpfr,
    phi_r: dict[int, mpfr],
    alpha_C: mpfr,
    *,
    prec: int,
) -> mpfr:
    """
    Valid cluster integral upper (L-0073R) replacing invalid B_C * max(phi):

      contrib <= |alpha| * B_C + sum_r |phi_r - alpha| * B_r

    At fixed shell factors phi_r (Route A weights). Do NOT multiply the full
    residual by max_phi (that revives the invalid B_C * max(phi) scaling).
    """
    _ctx(prec)
    max_phi = mpfr(0, precision=prec)
    per_shell = mpfr(0, precision=prec)
    residual = mul_up(abs_up_mpfr(alpha_C, prec), B_C)
    for r in shells:
        rf = phi_r[r]
        max_phi = max(max_phi, rf)
        diff = abs_up_mpfr(sub_up(rf, alpha_C), prec)
        residual = add_up(residual, mul_up(diff, B_r[r]))
        per_shell = add_up(per_shell, mul_up(B_r[r], rf))
    cluster_legacy = mul_up(B_C, max_phi)
    cluster_residual = residual
    # Valid: min(per-shell sum, lemma residual); legacy B_C*max_phi kept for audit only
    pick = per_shell if per_shell <= cluster_residual else cluster_residual
    return pick


def cluster_residual_vs_legacy_audit(
    shells: tuple[int, ...],
    B_r: dict[int, mpfr],
    B_C: mpfr,
    phi_r: dict[int, mpfr],
    *,
    prec: int,
) -> dict:
    """Compare L-0073R residual alpha choices vs invalid max*BC."""
    opt = cluster_alpha_optimize_contrib(shells, B_r, B_C, phi_r, prec=prec)
    rows = opt["candidates"]
    legacy = mul_up(B_C, max(phi_r.values(), key=float))
    best = opt["best"]
    return {
        "shells": list(shells),
        "legacy_invalid_hi": str(legacy),
        "residual_candidates": rows,
        "best_residual": best,
        "residual_tighter_than_legacy": float(best["contrib_hi"]) <= float(legacy),
        "alpha_scan": "knots_at_each_shell_phi_plus_heuristic_midpoints",
    }


def abs_up_mpfr(x: mpfr, prec: int) -> mpfr:
    with _ctx(prec, up=True):
        return mpfr(abs(x))


def cluster_alpha_candidates(phi_r: dict[int, mpfr]) -> dict[str, mpfr]:
    shells = sorted(phi_r.keys())
    if not shells:
        return {}
    vals = [phi_r[r] for r in shells]
    mid = vals[len(vals) // 2]
    cands: dict[str, mpfr] = {
        "min": min(vals, key=float),
        "max": max(vals, key=float),
        "midpoint": mid,
        "central_shell": phi_r[shells[len(shells) // 2]],
    }
    for r in shells:
        cands[f"shell_{r}"] = phi_r[r]
    return cands


def cluster_alpha_optimize_contrib(
    shells: tuple[int, ...],
    B_r: dict[int, mpfr],
    B_C: mpfr,
    phi_r: dict[int, mpfr],
    *,
    prec: int,
) -> dict:
    """Rigorous 1D scan: residual is piecewise linear; min at knot phi_r or endpoints."""
    cands = cluster_alpha_candidates(phi_r)
    rows = []
    for name, alpha in cands.items():
        c = cluster_residual_contrib_hi(shells, B_r, B_C, phi_r, alpha, prec=prec)
        rows.append({"alpha": name, "contrib_hi": str(c)})
    best = min(rows, key=lambda r: float(r["contrib_hi"]))
    return {"candidates": rows, "best": best}


def cluster_route_d_certificate(
    cluster_manifest: dict,
    *,
    per_shell_manifest: dict | None = None,
    prec: int | None = None,
) -> dict:
    """
    Route D with cluster blocks:
      I <= C(T)*T + sum_C min( sum_{r in C} B_r * int_phi(r), B_C * max int_phi )
    """
    prec = prec or cluster_manifest.get("arithmetic_backend", {}).get("precision_bits", 128)
    _ctx(int(prec))
    nu = mpfr(FROZEN.nu, precision=prec)
    T = mpfr(FROZEN.T, precision=prec)

    if per_shell_manifest:
        C_T_hi = mpfr(
            per_shell_manifest.get("full_best_hi", per_shell_manifest["full_frobenius_hi"])[
                "C_term_hi"
            ],
            precision=prec,
        )
    else:
        C_T_hi = mpfr(0, precision=prec)
        for blk in cluster_manifest["blocks"]:
            C_T_hi = add_up(C_T_hi, mpfr(blk["B_cluster_hi"], precision=prec))

    per_by: dict[int, mpfr] = {}
    if per_shell_manifest:
        for blk in per_shell_manifest.get("blocks", []):
            per_by[int(blk["shell"])] = mpfr(blk["C_term_hi"], precision=prec)

    integral_hi = mul_up(C_T_hi, T)
    for blk in cluster_manifest["blocks"]:
        B_C = mpfr(blk["B_cluster_hi"], precision=prec)
        per_sum = mpfr(0, precision=prec)
        max_phi = mpfr(0, precision=prec)
        for r in blk["shells"]:
            rf = mpfr(r, precision=prec)
            iphi = _shell_correction_hi(rf, nu=nu, T=T, prec=prec)
            max_phi = max(max_phi, iphi)
            B_r = per_by.get(int(r), B_C)
            per_sum = add_up(per_sum, mul_up(B_r, iphi))
        cluster_upper = mul_up(B_C, max_phi)
        contrib = per_sum if per_sum <= cluster_upper else cluster_upper
        integral_hi = add_up(integral_hi, contrib)

    I_lo, I_hi = i_star_interval(prec=int(prec))
    two_sqrt2 = sqrt_up(mpfr(8, precision=int(prec)))
    with _ctx(int(prec), up=True):
        floor_hi = mpfr(FROZEN.stokes_floor)
    omega_hi = add_up(floor_hi, div_up(integral_hi, two_sqrt2))

    return {
        "method": "route_D_cluster",
        "integral_hi": str(integral_hi),
        "omega_T_hi": str(omega_hi),
        "I_star_lo": str(I_lo),
        "closes_strict": bool(integral_hi < I_lo),
        "evidence_level": "N5",
    }


def build_cluster_manifest(
    n: int = 24,
    *,
    radii: tuple[int, ...] | None = None,
    n_clusters: int = 12,
    partition_strategy: str = "equal",
    clusters: list[ClusterSpec] | None = None,
    prec: int = 128,
    workers: int = 4,
    per_shell_manifest: dict | None = None,
    bound_method: str = "frobenius",
) -> dict:
    radii = radii or all_dealias_radii(n)
    if clusters is None:
        clusters = adaptive_cluster_partition(
            radii, partition_strategy, n_clusters=n_clusters
        )
    stc = shell_to_cluster_map(clusters)

    row_by: dict[int, list] | None = None
    col_by: dict[int, dict] | None = None
    if bound_method in ("best", "one_inf_only"):
        fro_cluster, row_by, col_by, n_cluster, D = _accumulate_cluster_bounds_single_pass(
            n, radii, stc, prec=prec, workers=workers
        )
    else:
        fro_cluster, n_cluster, D = _accumulate_frobenius_clusters_single_pass(
            n, radii, stc, prec=prec, workers=workers
        )

    per_by_shell: dict[int, mpfr] = {}
    if per_shell_manifest:
        for blk in per_shell_manifest.get("blocks", []):
            per_by_shell[int(blk["shell"])] = mpfr(blk["C_term_hi"], precision=prec)

    nu = mpfr(FROZEN.nu, precision=prec)
    T = mpfr(FROZEN.T, precision=prec)
    blocks = []
    for spec in clusters:
        cid = spec.cluster_id
        fro_sq = fro_cluster[cid]
        n_up = n_cluster[cid]
        bnd_fro = _bound_from_fro_sq(
            fro_sq, n_up, n, D, output_shell=spec.shells[0], prec=prec
        )
        if bound_method in ("best", "one_inf_only") and row_by is not None and col_by is not None:
            bnd_one = _bound_from_one_inf_rowlist(
                row_by[cid], col_by[cid], n_up, n, D, output_shell=spec.shells[0], prec=prec
            )
            if bound_method == "one_inf_only":
                bnd = bnd_one
            else:
                with _ctx(prec, up=True):
                    pick = (
                        bnd_one
                        if mpfr(bnd_one.C_term_hi) < mpfr(bnd_fro.C_term_hi)
                        else bnd_fro
                    )
                bnd = pick
            cluster_one_hi = bnd_one.C_term_hi
        else:
            bnd = bnd_fro
            cluster_one_hi = None

        B_C = mpfr(bnd.C_term_hi, precision=prec)
        phi_r: dict[int, mpfr] = {}
        B_r_map: dict[int, mpfr] = {}
        per_shell_sum = mpfr(0, precision=prec)
        max_factor = mpfr(0, precision=prec)
        per_detail = []
        for r in spec.shells:
            rf = mpfr(r, precision=prec)
            f_r = _shell_factor_hi(rf, nu=nu, T=T, prec=prec)
            phi_r[r] = f_r
            max_factor = max(max_factor, f_r)
            B_r = per_by_shell.get(r, B_C)
            B_r_map[r] = B_r
            per_shell_sum = add_up(per_shell_sum, mul_up(B_r, f_r))
            per_detail.append({"shell": r, "factor_hi": str(f_r), "B_r_hi": str(B_r)})

        legacy_invalid = mul_up(B_C, max_factor)
        if bound_method == "one_inf_only":
            audit = cluster_residual_vs_legacy_audit(
                tuple(spec.shells), B_r_map, B_C, phi_r, prec=prec
            )
            contrib_hi = mpfr(audit["best_residual"]["contrib_hi"], precision=prec)
            residual_method = f"L-0073R_{audit['best_residual']['alpha']}"
        else:
            cluster_upper = legacy_invalid
            contrib_hi = per_shell_sum if per_shell_sum <= cluster_upper else cluster_upper
            residual_method = "legacy_min"
            audit = None

        block = {
            **bnd.as_dict(),
            "cluster_id": cid,
            "shells": list(spec.shells),
            "B_cluster_hi": bnd.C_term_hi,
            "B_cluster_frobenius_hi": bnd_fro.C_term_hi,
            "max_factor_hi": str(max_factor),
            "per_shell_integral_hi": str(per_shell_sum),
            "cluster_integral_hi": str(legacy_invalid),
            "contrib_hi": str(contrib_hi),
            "method": f"route_a_prime_{residual_method}_{bound_method}",
            "per_shell_detail": per_detail,
        }
        if audit is not None:
            block["cluster_residual_L0073R"] = audit
        if cluster_one_hi is not None:
            block["B_cluster_one_inf_hi"] = cluster_one_hi
        blocks.append(block)

    manifest = {
        "n": n,
        "D_full": D,
        "n_clusters": len(clusters),
        "n_shells": len(set(radii)),
        "bound_method": bound_method,
        "partition": partition_strategy,
        "partition_strategy": partition_strategy,
        "clusters": [c.as_dict() for c in clusters],
        "blocks": blocks,
        "arithmetic_backend": ArithmeticBackend(precision_bits=prec).as_dict(),
        "convention": (
            "Route A': per cluster C, I_C from L-0073R residual when bound_method=one_inf_only; "
            "else legacy min(per-shell, B_C*max f). "
            f"B_C from {bound_method} combined-output bound."
        ),
    }
    payload = dict(manifest)
    payload.pop("manifest_sha256", None)
    manifest["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, default=str).encode()
    ).hexdigest()
    return manifest


def cluster_integral_certificate(
    manifest: dict,
    *,
    prec: int | None = None,
) -> dict:
    prec = prec or manifest.get("arithmetic_backend", {}).get("precision_bits", 128)
    _ctx(int(prec))
    integral_hi = mpfr(0, precision=prec)
    per_cluster = []
    for blk in manifest["blocks"]:
        c = mpfr(blk["contrib_hi"], precision=prec)
        integral_hi = add_up(integral_hi, c)
        per_cluster.append(
            {
                "cluster_id": blk["cluster_id"],
                "shells": blk["shells"],
                "contrib_hi": blk["contrib_hi"],
                "per_shell_integral_hi": blk["per_shell_integral_hi"],
                "cluster_integral_hi": blk["cluster_integral_hi"],
            }
        )

    I_lo, I_hi = i_star_interval(prec=int(prec))
    two_sqrt2 = sqrt_up(mpfr(8, precision=int(prec)))
    with _ctx(int(prec), up=True):
        floor_hi = mpfr(FROZEN.stokes_floor)
    omega_hi = add_up(floor_hi, div_up(integral_hi, two_sqrt2))
    with _ctx(int(prec), up=False):
        M_target_lo = mpfr(FROZEN.c0008_target_M)

    full_dealias = manifest.get("D_full") == FROZEN.D_full and manifest.get(
        "n_shells"
    ) == len(all_dealias_radii(manifest.get("n", 24)))
    closes = integral_hi < I_lo

    return {
        "method": "cluster_integral_route_A_prime",
        "n": manifest.get("n", 24),
        "evidence_level": "N5",
        "integral_hi": str(integral_hi),
        "I_star_lo": str(I_lo),
        "closes_strict": bool(closes),
        "omega_T_hi": str(omega_hi),
        "M_target_lo": str(M_target_lo),
        "closes_c0008": bool(closes and omega_hi < M_target_lo and full_dealias),
        "full_dealias": full_dealias,
        "n_clusters": manifest.get("n_clusters"),
        "per_cluster": per_cluster,
        "manifest_sha256": manifest.get("manifest_sha256"),
        "honesty": (
            "Per-cluster min(per-shell best sum, B_C max f(r)); B_C from "
            f"{manifest.get('bound_method', 'frobenius')} cluster combined bound. "
            "Triangle across clusters remains."
        ),
    }


def save_cluster_manifest(manifest: dict, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path
