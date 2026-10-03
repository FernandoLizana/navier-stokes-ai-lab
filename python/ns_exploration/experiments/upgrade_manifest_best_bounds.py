"""
Upgrade Frobenius shell manifest with per-shell 1-inf passes (memory-safe).

Resume layers:
  1. upgrade_best_checkpoint.json — completed shells (C_term one_inf)
  2. upgrade_shell_ck/shell_NNN.json — intra-shell pair progress (every 500k pairs)
  3. --watchdog — auto-restart on crash until all shells done
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from gmpy2 import mpfr

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.intervals import _ctx
from ns_exploration.terminal_weighted.opnorm_certified import (
    _SHELL_ONE_INF_CK_DIR,
    certified_L_op_one_inf_single_shell_hi,
)
from ns_exploration.terminal_weighted.shell_decomposition import save_shell_manifest


_CHECKPOINT = Path("experiments/terminal_weighted/upgrade_best_checkpoint.json")
_PARTIAL_OUT = Path("experiments/terminal_weighted/shell_manifest_best_partial.json")
_HEARTBEAT = Path("experiments/terminal_weighted/upgrade_best_heartbeat.txt")


def _touch_heartbeat(msg: str = "") -> None:
    _HEARTBEAT.parent.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}\n"
    _HEARTBEAT.write_text(line, encoding="utf-8")


def _worker_one_shell(args: tuple) -> dict:
    n, radii, shell, prec, ck_dir = args
    ck_path = Path(ck_dir) / f"shell_{shell:03d}.json"
    _touch_heartbeat(f"shell {shell} start")
    bnd = certified_L_op_one_inf_single_shell_hi(
        n,
        tuple(radii),
        shell,
        prec=prec,
        checkpoint_path=ck_path,
        heartbeat_path=_HEARTBEAT,
    )
    return {"shell": shell, **bnd.as_dict()}


def _load_done(manifest: dict, *, resume: bool) -> dict[int, dict]:
    if not resume or not _CHECKPOINT.is_file():
        return {}
    ck = json.loads(_CHECKPOINT.read_text(encoding="utf-8"))
    if ck.get("manifest_sha256") != manifest.get("manifest_sha256"):
        return {}
    return {int(k): v for k, v in ck.get("one_inf_by_shell", {}).items()}


def _save_ck(manifest: dict, done: dict[int, dict]) -> None:
    _CHECKPOINT.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "manifest_sha256": manifest.get("manifest_sha256"),
        "one_inf_by_shell": {str(k): v for k, v in sorted(done.items())},
        "n_completed": len(done),
        "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    tmp = _CHECKPOINT.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(_CHECKPOINT)


def _assemble_manifest(manifest: dict, done: dict[int, dict], *, prec: int) -> dict:
    blocks = []
    for blk in manifest["blocks"]:
        sh = int(blk["shell"])
        fro_hi = blk.get("C_term_frobenius_hi") or blk["C_term_hi"]
        one = done.get(sh, {})
        one_hi = one.get("C_term_hi", fro_hi)
        with _ctx(prec, up=True):
            pick_fro = mpfr(fro_hi)
            pick_one = mpfr(one_hi)
        if pick_one < pick_fro:
            chosen_hi, method = one_hi, "one_inf_mpfr_up"
        else:
            chosen_hi, method = fro_hi, "frobenius_mpfr_up"
        blocks.append(
            {
                **blk,
                "C_term_hi": chosen_hi,
                "C_term_frobenius_hi": fro_hi,
                "C_term_one_inf_hi": one_hi,
                "method": method,
                "n_updates_one_inf": one.get("n_updates", 0),
            }
        )
    upgraded = dict(manifest)
    upgraded["bound_method"] = "best"
    upgraded["blocks"] = blocks
    upgraded["full_best_hi"] = manifest.get("full_frobenius_hi") or manifest.get("full_best_hi")
    upgraded["upgrade_one_inf"] = {
        "shells_computed": len(done),
        "n_shells_total": manifest.get("n_shells"),
        "checkpoint": str(_CHECKPOINT),
        "shell_ck_dir": str(_SHELL_ONE_INF_CK_DIR),
        "method": "per_shell_one_inf_pass",
    }
    payload = dict(upgraded)
    payload.pop("manifest_sha256", None)
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()
    upgraded["manifest_sha256"] = digest
    return upgraded


def _save_partial(manifest: dict, done: dict[int, dict], *, prec: int) -> None:
    partial = _assemble_manifest(manifest, done, prec=prec)
    save_shell_manifest(partial, _PARTIAL_OUT)


def upgrade_manifest(
    manifest_path: str | Path,
    *,
    out_path: str | Path | None = None,
    prec: int = 128,
    workers: int = 1,
    resume: bool = True,
    shell_limit: int | None = None,
) -> dict:
    manifest_path = Path(manifest_path)
    out_path = Path(out_path or "experiments/terminal_weighted/shell_manifest_best_full.json")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    n = manifest["n"]
    radii = tuple(all_dealias_radii(n))
    done = _load_done(manifest, resume=resume)

    pending = [int(b["shell"]) for b in manifest["blocks"] if int(b["shell"]) not in done]
    if shell_limit is not None:
        pending = pending[:shell_limit]

    n_total = manifest.get("n_shells", len(manifest["blocks"]))
    print(
        f"upgrade best: {len(done)}/{n_total} done, {len(pending)} pending, workers={workers}",
        flush=True,
    )
    _touch_heartbeat(f"batch start done={len(done)} pending={len(pending)}")

    if not pending:
        upgraded = _assemble_manifest(manifest, done, prec=prec)
        save_shell_manifest(upgraded, out_path)
        return upgraded

    _SHELL_ONE_INF_CK_DIR.mkdir(parents=True, exist_ok=True)
    ck_dir = str(_SHELL_ONE_INF_CK_DIR)

    def _finish_shell(sh: int, res: dict) -> None:
        done[sh] = res
        _save_ck(manifest, done)
        _save_partial(manifest, done, prec=prec)
        _touch_heartbeat(f"shell {sh} done [{len(done)}/{n_total}]")
        print(f"  shell {sh} one_inf C_term_hi={res['C_term_hi']} [{len(done)}/{n_total}]", flush=True)

    if workers > 1:
        batch_size = workers
        for i in range(0, len(pending), batch_size):
            batch = pending[i : i + batch_size]
            tasks = [(n, list(radii), sh, prec, ck_dir) for sh in batch]
            with ProcessPoolExecutor(max_workers=len(batch)) as ex:
                futs = {ex.submit(_worker_one_shell, t): t[2] for t in tasks}
                for fut in as_completed(futs):
                    sh = futs[fut]
                    _finish_shell(sh, fut.result())
    else:
        for sh in pending:
            print(f"  computing 1-inf shell {sh}...", flush=True)
            ck_path = _SHELL_ONE_INF_CK_DIR / f"shell_{sh:03d}.json"
            bnd = certified_L_op_one_inf_single_shell_hi(
                n, radii, sh, prec=prec, checkpoint_path=ck_path, heartbeat_path=_HEARTBEAT
            )
            _finish_shell(sh, bnd.as_dict())

    upgraded = _assemble_manifest(manifest, done, prec=prec)
    save_shell_manifest(upgraded, out_path)
    print(f"saved {out_path} sha256={upgraded['manifest_sha256'][:16]}...", flush=True)
    return upgraded


def _count_completed(manifest_path: Path) -> tuple[int, int]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    done = _load_done(manifest, resume=True)
    total = manifest.get("n_shells", len(manifest["blocks"]))
    return len(done), total


def run_watchdog(
    manifest_path: str,
    *,
    out_path: str,
    prec: int,
    workers: int,
    retry_delay: float = 5.0,
    max_retries: int = 5000,
) -> None:
    mp = Path(manifest_path)
    for attempt in range(1, max_retries + 1):
        n_done, n_total = _count_completed(mp)
        if n_done >= n_total:
            print(f"watchdog: all {n_total} shells complete.", flush=True)
            upgrade_manifest(
                mp, out_path=out_path, prec=prec, workers=workers, resume=True
            )
            return
        print(
            f"watchdog: attempt {attempt} ({n_done}/{n_total} shells in checkpoint)...",
            flush=True,
        )
        try:
            upgrade_manifest(
                mp, out_path=out_path, prec=prec, workers=workers, resume=True
            )
            n_done, n_total = _count_completed(mp)
            if n_done >= n_total:
                print("watchdog: finished successfully.", flush=True)
                return
        except KeyboardInterrupt:
            print("watchdog: interrupted by user.", flush=True)
            raise
        except Exception:
            print(f"watchdog: crash on attempt {attempt}:", flush=True)
            traceback.print_exc()
        print(f"watchdog: sleeping {retry_delay}s before retry...", flush=True)
        time.sleep(retry_delay)
    raise RuntimeError(f"watchdog: exceeded {max_retries} retries")


def main(argv: list[str] | None = None) -> None:
    p = argparse.ArgumentParser(description="Upgrade Frobenius manifest with per-shell 1-inf")
    p.add_argument(
        "manifest",
        nargs="?",
        default="experiments/terminal_weighted/shell_manifest_frobenius_full.json",
    )
    p.add_argument("--out", default="experiments/terminal_weighted/shell_manifest_best_full.json")
    p.add_argument("--prec", type=int, default=128)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--no-resume", action="store_true")
    p.add_argument("--limit", type=int, default=None)
    p.add_argument(
        "--watchdog",
        action="store_true",
        help="auto-restart on crash until all shells complete",
    )
    p.add_argument("--retry-delay", type=float, default=5.0)
    p.add_argument(
        "--status",
        action="store_true",
        help="print completed/total shell count and exit",
    )
    args = p.parse_args(argv)

    if args.status:
        n_done, n_total = _count_completed(Path(args.manifest))
        print(f"{n_done}/{n_total}")
        return

    if args.watchdog:
        run_watchdog(
            args.manifest,
            out_path=args.out,
            prec=args.prec,
            workers=args.workers,
            retry_delay=args.retry_delay,
        )
    else:
        upgrade_manifest(
            args.manifest,
            out_path=args.out,
            prec=args.prec,
            workers=args.workers,
            resume=not args.no_resume,
            shell_limit=args.limit,
        )


if __name__ == "__main__":
    main()
