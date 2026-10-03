"""CLI: regenerate full-dealias manifest with one_inf_only per-shell passes."""

from __future__ import annotations

import argparse

from ns_exploration.terminal_weighted.repair_full_manifest import (
    build_repair_full_one_inf_manifest,
    repair_progress,
    run_watchdog,
)


def main() -> dict:
    p = argparse.ArgumentParser(description="C-0008 repair full-dealias one_inf manifest")
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--workers", type=int, default=1, help="parallel shells (max with --adaptive)")
    p.add_argument("--min-workers", type=int, default=1)
    p.add_argument("--max-workers", type=int, default=None)
    p.add_argument("--adaptive", action="store_true", help="RAM/OOM-aware worker cap + backoff on crash")
    p.add_argument("--shell-limit", type=int, default=None)
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--out", type=str, default=None)
    p.add_argument("--rational", action="store_true", help="mpfr_leray + mpfr_vdot_real coefficients")
    p.add_argument("--watchdog", action="store_true", help="auto-restart until all shells done")
    p.add_argument("--status", action="store_true", help="print progress and exit")
    p.add_argument("--retry-delay", type=float, default=10.0)
    args = p.parse_args()

    if args.status:
        prog = repair_progress(args.n, rational=args.rational)
        print(prog)
        return prog

    if args.watchdog:
        max_w = args.max_workers if args.max_workers is not None else args.workers
        manifest = run_watchdog(
            args.n,
            prec=args.prec,
            workers=max_w,
            retry_delay=args.retry_delay,
            adaptive=args.adaptive,
            min_workers=args.min_workers,
            max_workers=max_w,
            rational_coefficients=args.rational,
        )
    else:
        from pathlib import Path

        out = Path(args.out) if args.out else None
        manifest = build_repair_full_one_inf_manifest(
            args.n,
            prec=args.prec,
            workers=args.workers,
            resume=not args.no_resume,
            shell_limit=args.shell_limit,
            out_path=out,
            rational_coefficients=args.rational,
        )
    print(
        {
            "n": manifest["n"],
            "n_shells": manifest["n_shells"],
            "blocks": len(manifest["blocks"]),
            "complete": manifest.get("repair_complete"),
            "rational_coefficients": manifest.get("rational_coefficients", False),
            "sha256": manifest["manifest_sha256"],
        }
    )
    return manifest


if __name__ == "__main__":
    main()
