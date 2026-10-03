"""Watchdog: build L-0073R cluster manifest then cluster certificate."""

from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

PY = Path(__file__).resolve().parents[1]
REPO = PY.parent


def _cluster_path(n: int) -> Path:
    p = REPO / f"experiments/terminal_weighted/cluster_manifest_repair_full_n{n}_one_inf.json"
    if p.is_file():
        return p
    return PY / f"experiments/terminal_weighted/cluster_manifest_repair_full_n{n}_one_inf.json"


def main() -> dict:
    import argparse

    p = argparse.ArgumentParser()
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--workers", type=int, default=6)
    p.add_argument("--retry-delay", type=float, default=30.0)
    p.add_argument("--max-retries", type=int, default=20)
    args = p.parse_args()

    cluster = _cluster_path(args.n)
    for attempt in range(1, args.max_retries + 1):
        cluster = _cluster_path(args.n)
        if cluster.is_file():
            try:
                data = json.loads(cluster.read_text(encoding="utf-8"))
                if data.get("n_clusters", 0) >= 1 and len(data.get("blocks", [])) >= 1:
                    print(f"cluster manifest ready: {cluster}", flush=True)
                    break
            except Exception:
                pass
        print(f"watchdog cluster: attempt {attempt} building manifest...", flush=True)
        rc = subprocess.call(
            [
                sys.executable,
                str(PY / "scripts/build_c0008_repair_cluster_manifest.py"),
                "--n",
                str(args.n),
                "--workers",
                str(args.workers),
            ],
            cwd=str(PY),
        )
        if rc != 0:
            print(f"cluster build failed rc={rc}, retry in {args.retry_delay}s", flush=True)
            time.sleep(args.retry_delay)
            continue
        cluster = _cluster_path(args.n)
        if cluster.is_file():
            break
        time.sleep(args.retry_delay)
    else:
        raise RuntimeError("cluster manifest build exceeded retries")

    rc = subprocess.call(
        [
            sys.executable,
            str(PY / "scripts/build_c0008_repair_cluster_certificate.py"),
            "--n",
            str(args.n),
        ],
        cwd=str(PY),
    )
    if rc != 0:
        raise RuntimeError("cluster certificate build/verify failed")
    cert = REPO / "certificates/CERT-L0073R-repair-C0008-cluster-one-inf.json"
    rep = {"cluster_manifest": str(cluster), "cert": str(cert), "ok": cert.is_file()}
    print(json.dumps(rep, indent=2))
    return rep


if __name__ == "__main__":
    main()
