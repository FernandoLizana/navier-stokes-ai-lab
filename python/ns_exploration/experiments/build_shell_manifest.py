"""Build shell manifest only (Frobenius pass + save)."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.shell_decomposition import (
    build_shell_manifest,
    save_shell_manifest,
)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=("band", "full"), default="full")
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--bound-method", choices=("frobenius", "one_inf", "best"), default="best")
    args = p.parse_args()
    radii = tuple(range(1, 7)) if args.mode == "band" else all_dealias_radii(24)
    print(f"manifest mode={args.mode} radii={len(radii)} workers={args.workers}", flush=True)
    manifest = build_shell_manifest(24, prec=args.prec, radii=radii, workers=args.workers, bound_method=args.bound_method)
    path = Path("experiments/terminal_weighted/shell_manifest.json")
    save_shell_manifest(manifest, path)
    print(f"saved {path} D={manifest['D_full']} shells={manifest['n_shells']}", flush=True)
    print(f"sha256={manifest.get('manifest_sha256')}", flush=True)


if __name__ == "__main__":
    main()
