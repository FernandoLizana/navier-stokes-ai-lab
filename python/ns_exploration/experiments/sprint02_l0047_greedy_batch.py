"""
Sprint helper: batch-greedy enlarge SUPPORT_CR0011 under float32 dense-M cap.

Exploration only (N2/N7 hybrid). Package L-0047 only if an enlarged support
satisfies C_fullsym ≤ C_† with documented dtype/margin.

Writes an incremental checkpoint after every accept/reject so a crash/hibernate
does not lose progress. Resume with --resume <checkpoint.json>.

FINITE Galerkin only. Not continuum. Not Clay.
"""

from __future__ import annotations

import argparse
import json
import os
import time

# Cap BLAS/OpenMP threads BEFORE importing numpy so a background exploration
# does not saturate every core (this is what froze the machine). Override with
# NS_THREADS=<k>. Default: leave 2 cores free.
_default_threads = max(1, (os.cpu_count() or 4) - 2)
_ns_threads = os.environ.get("NS_THREADS", str(_default_threads))
for _var in (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
):
    os.environ.setdefault(_var, _ns_threads)

from pathlib import Path

import numpy as np


def _lower_priority() -> None:
    """Best-effort: drop process priority so the UI stays responsive."""
    try:
        import psutil  # type: ignore

        p = psutil.Process()
        if os.name == "nt":
            p.nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
        else:
            p.nice(10)
    except Exception:
        try:
            os.nice(10)  # POSIX fallback
        except (AttributeError, OSError):
            pass


from ns_exploration.conjectures.l0021_quartic_shell import shell_modes_by_r
from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0045_streaming_fullsym import streaming_C_fullsym
from ns_exploration.conjectures.l0046_greedy_fullsym import SUPPORT_CR0011
from ns_exploration.experiments.sprintB_physical_tensor import build_hermitian_basis


def _shell_D(n: int, r: int) -> int:
    return build_hermitian_basis(n, (r,)).D


def _mem_gb(D: int, itemsize: int = 4) -> float:
    return D * (D * (D + 1) // 2) * itemsize / 1e9


def _save_ckpt(path: Path | None, payload: dict) -> None:
    if path is None:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    payload = {**payload, "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)


def run_batch_greedy(
    n: int = 24,
    mem_gb: float = 12.0,
    batch: int = 6,
    seed_extra: tuple[int, ...] = (),
    skip_rejected: tuple[int, ...] = (),
    tol: float = 1e-5,
    ckpt_path: Path | None = None,
) -> dict:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    by = shell_modes_by_r(n)
    support = sorted(set(SUPPORT_CR0011) | set(seed_extra))
    rejected: list[tuple[int, float]] = []
    skip = set(skip_rejected)
    missing = sorted(
        (r for r in by if r not in set(support) and r not in skip),
        key=lambda r: (_shell_D(n, r), r),
    )
    print(
        f"Cd={Cd} start n={len(support)} D={build_hermitian_basis(n, tuple(support)).D} "
        f"missing={len(missing)} skip={sorted(skip)}",
        flush=True,
    )
    s0 = streaming_C_fullsym(n, tuple(support), dtype=np.float32)
    Ccur = float(s0["C_fullsym"])
    print(f"start C={Ccur}", flush=True)

    added: list[int] = []

    def snapshot(done: bool = False) -> dict:
        return {
            "done": done,
            "C_dagger": Cd,
            "C_fullsym": Ccur,
            "support": support,
            "added": added,
            "rejected": rejected,
            "skip_rejected": sorted(skip),
            "n_shells": len(support),
            "D": build_hermitian_basis(n, tuple(support)).D,
            "missing_left": [r for r in sorted(by) if r not in set(support)],
            "closes_c0007_all_ic": all(r in set(support) for r in by),
            "dtype": "float32",
            "mem_gb_cap": mem_gb,
        }

    _save_ckpt(ckpt_path, snapshot(False))

    i = 0
    while i < len(missing):
        cand: list[int] = []
        while i < len(missing) and len(cand) < batch:
            r = missing[i]
            trial = tuple(sorted(support + cand + [r]))
            D = build_hermitian_basis(n, trial).D
            if _mem_gb(D) > mem_gb:
                if not cand:
                    print(f"skip r={r} D={D} M32_GB={_mem_gb(D):.1f}", flush=True)
                    i += 1
                break
            cand.append(r)
            i += 1
        if not cand:
            if i >= len(missing):
                break
            continue
        trial = tuple(sorted(support + cand))
        D = build_hermitian_basis(n, trial).D
        print(f"try batch={cand} D={D} M32_GB={_mem_gb(D):.2f}...", flush=True)
        s = streaming_C_fullsym(n, trial, dtype=np.float32)
        C = float(s["C_fullsym"])
        ok = C <= Cd + tol
        print(f"  C={C:.6f} ok={ok}", flush=True)
        if ok:
            support = sorted(support + cand)
            Ccur = C
            added.extend(cand)
            print(f"  ADDED n={len(support)} C={Ccur:.6f}", flush=True)
            _save_ckpt(ckpt_path, snapshot(False))
            continue
        print("  batch fail -> singles", flush=True)
        for r in cand:
            trial = tuple(sorted(support + [r]))
            D = build_hermitian_basis(n, trial).D
            if _mem_gb(D) > mem_gb:
                print(f"  skip r={r}", flush=True)
                continue
            print(f"  try r={r} D={D}...", flush=True)
            s = streaming_C_fullsym(n, trial, dtype=np.float32)
            C = float(s["C_fullsym"])
            ok = C <= Cd + tol
            print(f"    C={C:.6f} ok={ok}", flush=True)
            if ok:
                support = sorted(support + [r])
                Ccur = C
                added.append(r)
                print(f"    ADDED n={len(support)}", flush=True)
            else:
                rejected.append((r, C))
                skip.add(r)
                print("    reject", flush=True)
            _save_ckpt(ckpt_path, snapshot(False))

    out = snapshot(True)
    print("FINAL", json.dumps({k: out[k] for k in out if k != "support"}), flush=True)
    print("FINAL support", support, flush=True)
    _save_ckpt(ckpt_path, out)
    return out


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--mem-gb", type=float, default=12.0)
    p.add_argument("--batch", type=int, default=6)
    p.add_argument("--seed-extra", type=int, nargs="*", default=[43, 44])
    p.add_argument(
        "--skip-rejected",
        type=int,
        nargs="*",
        default=[],
        help="Shells already known to exceed C_dagger (skip recompute).",
    )
    p.add_argument(
        "--resume",
        type=str,
        default="",
        help="Load support/rejected from a previous checkpoint JSON.",
    )
    p.add_argument(
        "--out",
        type=str,
        default="reports/l0047_greedy_batch_result.json",
    )
    args = p.parse_args()
    _lower_priority()
    print(f"threads={_ns_threads} cpu_count={os.cpu_count()}", flush=True)

    seed = list(args.seed_extra)
    skip = list(args.skip_rejected)
    if args.resume:
        prev = json.loads(Path(args.resume).read_text(encoding="utf-8"))
        # seed_extra = support minus CR0011 baseline
        seed = [r for r in prev.get("support", []) if r not in SUPPORT_CR0011]
        skip = sorted(
            set(skip)
            | {int(x[0]) if isinstance(x, (list, tuple)) else int(x) for x in prev.get("rejected", [])}
            | set(prev.get("skip_rejected", []))
        )
        print(f"resume from {args.resume}: n={len(prev.get('support', []))} skip={skip}", flush=True)

    out_path = Path(args.out)
    out = run_batch_greedy(
        mem_gb=args.mem_gb,
        batch=args.batch,
        seed_extra=tuple(seed),
        skip_rejected=tuple(skip),
        ckpt_path=out_path,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(f"wrote {out_path}", flush=True)


if __name__ == "__main__":
    main()
