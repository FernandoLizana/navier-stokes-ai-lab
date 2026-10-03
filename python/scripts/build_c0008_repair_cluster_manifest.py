"""Build L-0073R cluster manifest from repair per-shell manifest."""

from __future__ import annotations

import json
from pathlib import Path

from ns_exploration.terminal_weighted.cluster_bounds import save_cluster_manifest
from ns_exploration.terminal_weighted.portable_paths import data_root, rational_tag
from ns_exploration.validation.c0008_repair_pipeline import build_cluster_manifest_repair

PY = Path(__file__).resolve().parents[1]
REPO = PY.parent


def main() -> dict:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--n-clusters", type=int, default=12)
    p.add_argument("--rational", action="store_true")
    p.add_argument("--manifest", type=str, default=None)
    p.add_argument("--out", type=str, default=None)
    args = p.parse_args()

    tag = rational_tag(args.rational)
    manifest_path = (
        Path(args.manifest)
        if args.manifest
        else data_root() / f"shell_manifest_repair_full_one_inf_n{args.n}{tag}.json"
    )
    if not manifest_path.is_file():
        manifest_path = PY / f"experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n{args.n}{tag}.json"
    out = (
        Path(args.out)
        if args.out
        else data_root() / f"cluster_manifest_repair_full_n{args.n}_one_inf{tag}.json"
    )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    cluster = build_cluster_manifest_repair(
        manifest, prec=args.prec, workers=args.workers, n_clusters=args.n_clusters
    )
    save_cluster_manifest(cluster, out)
    rep = {
        "manifest": str(manifest_path),
        "cluster_manifest": str(out),
        "n_clusters": cluster["n_clusters"],
        "sha256": cluster["manifest_sha256"],
    }
    print(json.dumps(rep, indent=2))
    return rep


if __name__ == "__main__":
    main()
