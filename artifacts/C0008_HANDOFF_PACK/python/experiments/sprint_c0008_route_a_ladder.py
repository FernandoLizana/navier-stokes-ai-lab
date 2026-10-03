"""Route A ladder: adaptive cluster partitions + Route A'/D sweep."""

from __future__ import annotations

import argparse
import json
import sys
import time
import traceback
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.cluster_bounds import (
    build_cluster_manifest,
    cluster_integral_certificate,
    cluster_route_d_certificate,
    save_cluster_manifest,
)
from ns_exploration.terminal_weighted.constants import FROZEN
from ns_exploration.terminal_weighted.intervals import i_star_interval

ROOT = Path(__file__).resolve().parents[3]

FULL_SWEEP_SPECS: tuple[tuple[str, str, int], ...] = (
    ("equal_12", "equal", 12),
    ("equal_24", "equal", 24),
    ("equal_48", "equal", 48),
    ("fine_low", "fine_low", 12),
)

L0073_MANIFEST = ROOT / "experiments/terminal_weighted/shell_manifest_cluster_best_full.json"


def _best_per_shell_path(mode: str) -> Path:
    if mode == "band":
        return ROOT / "experiments/terminal_weighted/shell_manifest.json"
    return ROOT / "experiments/terminal_weighted/shell_manifest_best_full.json"


def _checkpoint_path(mode: str) -> Path:
    return ROOT / f"experiments/terminal_weighted/route_a_ladder_{mode}_checkpoint.json"


def _manifest_cache_path(label: str) -> Path:
    safe = label.replace("/", "_")
    return ROOT / f"experiments/terminal_weighted/route_a_manifest_{safe}.json"


def _row_from_manifest(
    label: str,
    strategy: str,
    manifest: dict,
    *,
    per_shell: dict,
    I_lo: mpfr,
    runtime_s: float,
) -> dict:
    route_a = cluster_integral_certificate(manifest, prec=manifest["arithmetic_backend"]["precision_bits"])
    route_d = cluster_route_d_certificate(
        manifest, per_shell_manifest=per_shell, prec=manifest["arithmetic_backend"]["precision_bits"]
    )
    I_a = mpfr(route_a["integral_hi"])
    I_d = mpfr(route_d["integral_hi"])
    pick = "A" if I_a <= I_d else "D"
    I_best = I_a if pick == "A" else I_d
    omega = float(route_a["omega_T_hi"] if pick == "A" else route_d["omega_T_hi"])
    return {
        "label": label,
        "strategy": strategy,
        "n_clusters": manifest["n_clusters"],
        "route_a_I_hi": str(I_a),
        "route_d_I_hi": str(I_d),
        "best_route": pick,
        "best_I_hi": str(I_best),
        "omega_T_hi": omega,
        "closes_c0008": float(I_best) < float(I_lo) and omega < FROZEN.c0008_target_M,
        "runtime_s": round(runtime_s, 2),
        "manifest_path": str(_manifest_cache_path(label).relative_to(ROOT)).replace("\\", "/"),
    }


def seed_equal_12_from_l0073(*, prec: int = 128) -> dict | None:
    """Import equal_12 row from precomputed L-0073 manifest (no D^2 pass)."""
    if not L0073_MANIFEST.is_file():
        return None
    per_shell = json.loads(_best_per_shell_path("full").read_text(encoding="utf-8"))
    manifest = json.loads(L0073_MANIFEST.read_text(encoding="utf-8"))
    I_lo, _ = i_star_interval(prec=prec)
    man_path = _manifest_cache_path("equal_12")
    save_cluster_manifest(manifest, man_path)
    row = _row_from_manifest("equal_12", "equal", manifest, per_shell=per_shell, I_lo=I_lo, runtime_s=0.0)
    row["seeded_from"] = "L-0073"
    return row


def _save_checkpoint(mode: str, *, rows: list[dict], best: dict | None) -> None:
    ck = {
        "mode": mode,
        "rows": rows,
        "best": {k: v for k, v in (best or {}).items() if k != "manifest"},
        "best_manifest_path": best.get("manifest_path") if best else None,
    }
    _checkpoint_path(mode).write_text(json.dumps(ck, indent=2), encoding="utf-8")


def _load_checkpoint(mode: str) -> dict | None:
    path = _checkpoint_path(mode)
    if not path.is_file():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def finalize_ladder_summary(
    *,
    mode: str,
    rows: list[dict],
    best: dict | None,
    prec: int = 128,
) -> dict:
    I_lo, _ = i_star_interval(prec=prec)
    if best is None and rows:
        best = min(rows, key=lambda r: float(r["best_I_hi"]))

    suffix = mode
    summary = {
        "mode": mode,
        "rows": rows,
        "best": best,
        "best_manifest": best.get("manifest_path") if best else None,
        "I_star_lo": str(I_lo),
        "M_target": FROZEN.c0008_target_M,
        "complete": len(rows) >= len(FULL_SWEEP_SPECS) if mode == "full" else True,
    }
    # complete flag may be overridden by caller after all specs succeed

    out = ROOT / f"experiments/terminal_weighted/route_a_ladder_{suffix}.json"
    out.write_text(json.dumps(summary, indent=2), encoding="utf-8")

    report = ROOT / f"reports/C0008_ROUTE_A_LADDER_{suffix.upper()}.md"
    report.write_text(_render(summary), encoding="utf-8")

    if best and best.get("manifest_path"):
        src = ROOT / best["manifest_path"]
        if src.is_file():
            dst = ROOT / f"experiments/terminal_weighted/shell_manifest_route_a_best_{suffix}.json"
            dst.write_text(src.read_text(encoding="utf-8"), encoding="utf-8")
            summary["best_manifest"] = str(dst.relative_to(ROOT)).replace("\\", "/")

    return summary


def run_ladder(
    *,
    mode: str = "band",
    strategies: tuple[str, ...] = ("equal", "fine_low", "singleton"),
    cluster_counts: tuple[int, ...] = (2, 4, 6, 11, 24),
    prec: int = 128,
    workers: int | None = None,
    resume: bool = False,
    only_labels: tuple[str, ...] | None = None,
    seed_l0073: bool = False,
) -> dict:
    if workers is None:
        workers = 1 if mode == "band" else 4
    n = 24
    radii = tuple(range(1, 7)) if mode == "band" else all_dealias_radii(n)
    per_shell = json.loads(_best_per_shell_path(mode).read_text(encoding="utf-8"))
    I_lo, _ = i_star_interval(prec=prec)

    rows: list[dict] = []
    best: dict | None = None
    done_labels: set[str] = set()

    if resume:
        ck = _load_checkpoint(mode)
        if ck:
            rows = [r for r in ck.get("rows", []) if "error" not in r]
            done_labels = {r["label"] for r in rows}
            if ck.get("best"):
                best = ck["best"]
            print(f"  resume: {len(rows)} configs already done", flush=True)

    if seed_l0073 and mode == "full" and "equal_12" not in done_labels:
        seeded = seed_equal_12_from_l0073(prec=prec)
        if seeded:
            rows.append(seeded)
            done_labels.add("equal_12")
            if best is None or mpfr(seeded["best_I_hi"]) < mpfr(best["best_I_hi"]):
                best = seeded
            _save_checkpoint(mode, rows=rows, best=best)
            print(
                f"  equal_12 (seeded L-0073): A={float(seeded['route_a_I_hi']):.4f} "
                f"D={float(seeded['route_d_I_hi']):.4f} omega={seeded['omega_T_hi']:.4f}",
                flush=True,
            )

    specs: list[tuple[str, str, int]] = []
    if mode == "full" and strategies == ("equal", "fine_low") and cluster_counts == (12, 24, 48):
        specs = list(FULL_SWEEP_SPECS)
    else:
        for strategy in strategies:
            counts = cluster_counts if strategy == "equal" else (12,)
            for nc in counts:
                label = f"{strategy}_{nc}" if strategy == "equal" else strategy
                specs.append((label, strategy, nc if strategy == "equal" else 12))

    if only_labels:
        specs = [s for s in specs if s[0] in only_labels]

    for label, strategy, nc in specs:
        if label in done_labels:
            continue
        t0 = time.perf_counter()
        print(f"  {label}...", flush=True)
        try:
            manifest = build_cluster_manifest(
                n,
                radii=radii,
                n_clusters=nc,
                partition_strategy=strategy,
                prec=prec,
                workers=workers,
                per_shell_manifest=per_shell,
                bound_method="best",
            )
        except Exception as exc:
            err = str(exc) or repr(exc)
            tb = traceback.format_exc()
            print(f"    ERROR {label}: {err}", flush=True)
            rows.append({"label": label, "error": err, "traceback": tb})
            _save_checkpoint(mode, rows=rows, best=best)
            continue

        man_path = _manifest_cache_path(label)
        save_cluster_manifest(manifest, man_path)
        row = _row_from_manifest(
            label, strategy, manifest, per_shell=per_shell, I_lo=I_lo, runtime_s=time.perf_counter() - t0
        )
        rows.append(row)
        if best is None or mpfr(row["best_I_hi"]) < mpfr(best["best_I_hi"]):
            best = row
        _save_checkpoint(mode, rows=rows, best=best)
        print(
            f"    A={float(row['route_a_I_hi']):.4f} D={float(row['route_d_I_hi']):.4f} "
            f"pick={row['best_route']} omega={row['omega_T_hi']:.4f}",
            flush=True,
        )

    if not rows:
        raise RuntimeError("no ladder rows produced")

    ok_rows = [r for r in rows if "error" not in r]
    if mode == "full":
        complete = all(any(r["label"] == lab for r in ok_rows) for lab, _, _ in FULL_SWEEP_SPECS)
    else:
        complete = len(ok_rows) == len(specs) if specs else True

    summary = finalize_ladder_summary(mode=mode, rows=ok_rows, best=best, prec=prec)
    summary["complete"] = complete
    (ROOT / f"experiments/terminal_weighted/route_a_ladder_{mode}.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    if _checkpoint_path(mode).is_file() and complete:
        _checkpoint_path(mode).unlink()
    if not complete and not only_labels:
        pending = [lab for lab, _, _ in FULL_SWEEP_SPECS if lab not in {r["label"] for r in ok_rows}]
        print(f"  incomplete: pending {pending}", flush=True)
        sys.exit(2)
    return summary


def _render(summary: dict) -> str:
    status = "complete" if summary.get("complete") else "partial"
    lines = [
        f"# C-0008 Route A Ladder ({summary['mode']}) — {status}",
        "",
    ]
    if summary.get("best"):
        b = summary["best"]
        lines.append(
            f"**Best:** `{b['label']}` — I_hi={float(b['best_I_hi']):.4f}, "
            f"Omega={b['omega_T_hi']:.4f}, route {b['best_route']}"
        )
        lines.append("")
    lines.extend(
        [
            "| Label | clusters | Route A | Route D | best | Omega |",
            "|-------|----------|---------|---------|------|-------|",
        ]
    )
    for r in summary["rows"]:
        if "error" in r:
            lines.append(f"| {r['label']} | - | ERR | ERR | - | - |")
            continue
        lines.append(
            f"| {r['label']} | {r['n_clusters']} | {float(r['route_a_I_hi']):.4f} | "
            f"{float(r['route_d_I_hi']):.4f} | {r['best_route']} | {r['omega_T_hi']:.4f} |"
        )
    return "\n".join(lines)


def main() -> None:
    p = argparse.ArgumentParser(description="Route A cluster ladder sweep")
    p.add_argument("--mode", choices=("band", "full"), default="band")
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--workers", type=int, default=None)
    p.add_argument("--resume", action="store_true", help="skip configs in checkpoint")
    p.add_argument(
        "--seed-l0073",
        action="store_true",
        help="import equal_12 from L-0073 manifest (no recompute)",
    )
    p.add_argument(
        "--only",
        nargs="+",
        metavar="LABEL",
        help="run only these labels (e.g. equal_24 equal_48 fine_low)",
    )
    p.add_argument(
        "--full-sweep",
        action="store_true",
        help="full dealias: equal 12,24,48 + fine_low",
    )
    args = p.parse_args()
    if args.workers is not None:
        workers = args.workers
    elif args.mode == "band" or (not args.full_sweep and args.mode != "full"):
        workers = 1
    elif sys.platform == "win32":
        workers = 2
    else:
        workers = 4
    if args.mode == "full" or args.full_sweep:
        strategies = ("equal", "fine_low")
        counts = (12, 24, 48)
    else:
        strategies = ("equal", "fine_low", "singleton")
        counts = (2, 4, 6, 11)
    s = run_ladder(
        mode=args.mode if not args.full_sweep else "full",
        strategies=strategies,
        cluster_counts=counts,
        prec=args.prec,
        workers=workers,
        resume=args.resume,
        only_labels=tuple(args.only) if args.only else None,
        seed_l0073=args.seed_l0073 and not args.only,
    )
    print(json.dumps(s, indent=2))


if __name__ == "__main__":
    main()
