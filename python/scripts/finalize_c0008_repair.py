"""Finalize C-0008 repair: enrich manifest, cluster L-0073R, certificate, verify."""

from __future__ import annotations

import argparse
from pathlib import Path

from ns_exploration.validation.c0008_repair_pipeline import run_finalize_pipeline


def main() -> dict:
    p = argparse.ArgumentParser(description="C-0008 repair finalize pipeline")
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--workers", type=int, default=1)
    p.add_argument("--manifest", type=str, default=None)
    p.add_argument("--band", action="store_true", help="use band6 repair manifest")
    args = p.parse_args()

    repo = Path(__file__).resolve().parents[2]
    py_exp = Path(__file__).resolve().parents[1] / "experiments/terminal_weighted"
    if args.band:
        manifest_path = repo / "experiments/terminal_weighted/shell_manifest_repair_band6_one_inf.json"
        out_cert = repo / "certificates/CERT-C0008-repair-band6-one-inf.json"
        out_cluster = repo / "experiments/terminal_weighted/cluster_manifest_repair_band6_one_inf.json"
    else:
        full_repo = repo / f"experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n{args.n}.json"
        full_py = py_exp / f"shell_manifest_repair_full_one_inf_n{args.n}.json"
        partial_py = py_exp / f"shell_manifest_repair_full_one_inf_n{args.n}_partial.json"
        if args.manifest:
            manifest_path = Path(args.manifest)
        elif full_py.is_file():
            manifest_path = full_py
        elif full_repo.is_file():
            manifest_path = full_repo
        else:
            manifest_path = partial_py if partial_py.is_file() else full_repo
        out_cert = repo / f"certificates/CERT-C0008-repair-full-n{args.n}-one-inf.json"
        out_cluster = repo / f"experiments/terminal_weighted/cluster_manifest_repair_full_n{args.n}_one_inf.json"

    rep = run_finalize_pipeline(
        manifest_path,
        prec=args.prec,
        workers=args.workers,
        out_cert=out_cert,
        out_cluster=out_cluster,
    )
    print(rep)
    return rep


if __name__ == "__main__":
    main()
