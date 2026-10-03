"""
Build full-dealias shell manifest with one_inf_only per-shell passes (repair baseline).

Checkpoints (per n): repair_one_inf_checkpoint_n{N}.json
Partial output: shell_manifest_repair_full_one_inf_n{N}_partial.json
Intra-shell ck: shell_one_inf_ck/n{N}/shell_NNN.json
Watchdog state: repair_watchdog_state_n{N}.json
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
import time
import traceback
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

from ns_exploration.conjectures.l0048_matrixfree_fullsym import all_dealias_radii
from ns_exploration.terminal_weighted.intervals import ArithmeticBackend
from ns_exploration.terminal_weighted.portable_paths import data_root, rational_tag
from ns_exploration.terminal_weighted.opnorm_certified import (
    certified_L_op_one_inf_single_shell_hi,
)

# Conservative RAM budget per parallel shell worker (MPFR + basis + dict accumulators).
_RAM_GB_PER_WORKER = 1.75
_OOM_MARKERS = (
    "memoryerror",
    "out of memory",
    "cannot allocate",
    "bad_alloc",
    "terminating process due to low memory",
    "page file",
)


def _available_ram_gb() -> float:
    """Best-effort available physical RAM in GB."""
    if sys.platform == "win32":
        import ctypes

        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        stat = MEMORYSTATUSEX()
        stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        if ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat)):
            return stat.ullAvailPhys / (1024**3)
    try:
        import psutil  # type: ignore

        return psutil.virtual_memory().available / (1024**3)
    except Exception:
        pass
    return 8.0


def _is_memory_error(exc: BaseException) -> bool:
    if isinstance(exc, MemoryError):
        return True
    msg = str(exc).lower()
    return any(m in msg for m in _OOM_MARKERS)


def _rational_tag(rational: bool) -> str:
    return rational_tag(rational)


def _watchdog_state_path(n: int, *, rational: bool = False) -> Path:
    tag = _rational_tag(rational)
    return data_root() / f"repair_watchdog_state_n{n}{tag}.json"


def _load_watchdog_state(n: int, *, rational: bool = False) -> dict:
    path = _watchdog_state_path(n, rational=rational)
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}


def _save_watchdog_state(n: int, state: dict, *, rational: bool = False) -> None:
    path = _watchdog_state_path(n, rational=rational)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(state, indent=2), encoding="utf-8")
    tmp.replace(path)


def pick_workers(
    requested: int,
    *,
    n_pending: int = 1,
    ram_gb_per_worker: float = _RAM_GB_PER_WORKER,
) -> int:
    """Cap parallel shells by RAM, CPU count, and pending work."""
    cpu = os.cpu_count() or 4
    # Leave ~2 GB + one core for OS / parent process.
    ram_cap = max(1, int((_available_ram_gb() - 2.0) // ram_gb_per_worker))
    cpu_cap = max(1, cpu - 1)
    cap = min(requested, ram_cap, cpu_cap, max(1, n_pending))
    return max(1, cap)


def _backoff_workers(current: int, exc: BaseException | None, *, min_workers: int) -> int:
    if exc is None:
        return current
    if _is_memory_error(exc):
        return max(min_workers, current // 2)
    # Generic crash: step down one notch (thermal / process pool instability).
    return max(min_workers, current - 1)


def _paths(n: int, *, rational: bool = False) -> tuple[Path, Path, Path, Path]:
    tag = _rational_tag(rational)
    base = data_root()
    return (
        base / f"repair_one_inf_checkpoint_n{n}{tag}.json",
        base / f"shell_manifest_repair_full_one_inf_n{n}{tag}_partial.json",
        base / f"repair_one_inf_heartbeat_n{n}{tag}.txt",
        base / "shell_one_inf_ck" / f"n{n}{tag}",
    )


def _resolve_shell_ck_path(n: int, shell: int, ck_shell_dir: Path, *, rational: bool = False) -> Path:
    """Use n-specific ck; migrate legacy upgrade_shell_ck if present (float path only)."""
    ck_shell_dir.mkdir(parents=True, exist_ok=True)
    target = ck_shell_dir / f"shell_{shell:03d}.json"
    if target.is_file():
        return target
    if rational:
        return target
    legacy = data_root() / "upgrade_shell_ck" / f"shell_{shell:03d}.json"
    if not legacy.is_file():
        legacy = Path("experiments/terminal_weighted/upgrade_shell_ck") / f"shell_{shell:03d}.json"
    if legacy.is_file():
        ck = json.loads(legacy.read_text(encoding="utf-8"))
        if ck.get("n") == n and int(ck.get("shell", -1)) == shell:
            target.write_text(legacy.read_text(encoding="utf-8"), encoding="utf-8")
            print(f"    migrated shell {shell} ck from {legacy} ({ck.get('done_pairs')}/{ck.get('n_pairs')})", flush=True)
    return target


def _touch_heartbeat(path: Path, msg: str = "") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = f"{time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())} {msg}\n"
    path.write_text(line, encoding="utf-8")


def _worker_one_shell(args: tuple) -> dict:
    n, radii, shell, prec, ck_dir, hb_path, rational = args
    ck_path = _resolve_shell_ck_path(n, shell, Path(ck_dir), rational=rational)
    _touch_heartbeat(Path(hb_path), f"shell {shell} start n={n} rational={rational}")
    bnd = certified_L_op_one_inf_single_shell_hi(
        n,
        tuple(radii),
        shell,
        prec=prec,
        checkpoint_path=ck_path,
        heartbeat_path=Path(hb_path),
        rational_coefficients=rational,
    )
    return {"shell": shell, **bnd.as_dict(), "C_term_one_inf_hi": bnd.C_term_hi}


def _load_done(
    *, n: int, radii: tuple[int, ...], prec: int, resume: bool, rational: bool = False
) -> dict[int, dict]:
    ck_path, _, _, _ = _paths(n, rational=rational)
    legacy = data_root() / "repair_one_inf_checkpoint.json"
    if not legacy.is_file():
        legacy = Path("experiments/terminal_weighted/repair_one_inf_checkpoint.json")
    if resume and not ck_path.is_file() and legacy.is_file():
        old = json.loads(legacy.read_text(encoding="utf-8"))
        if old.get("n") == n and tuple(old.get("radii", [])) == radii and old.get("prec") == prec:
            ck_path.parent.mkdir(parents=True, exist_ok=True)
            ck_path.write_text(json.dumps(old, indent=2), encoding="utf-8")
    if not resume or not ck_path.is_file():
        return {}
    ck = json.loads(ck_path.read_text(encoding="utf-8"))
    if ck.get("n") != n or tuple(ck.get("radii", [])) != radii or ck.get("prec") != prec:
        return {}
    if ck.get("rational_coefficients", False) != rational:
        return {}
    return {int(k): v for k, v in ck.get("blocks_by_shell", {}).items()}


def _save_ck(
    *, n: int, radii: tuple[int, ...], prec: int, done: dict[int, dict], rational: bool = False
) -> None:
    ck_path, _, _, _ = _paths(n, rational=rational)
    ck_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "n": n,
        "radii": list(radii),
        "prec": prec,
        "rational_coefficients": rational,
        "blocks_by_shell": {str(k): v for k, v in sorted(done.items())},
        "n_completed": len(done),
        "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    tmp = ck_path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(ck_path)


def assemble_repair_manifest(
    *,
    n: int,
    radii: tuple[int, ...],
    blocks_by_shell: dict[int, dict],
    prec: int,
    complete: bool,
    rational: bool = False,
) -> dict:
    shells = sorted(set(radii))
    blocks = []
    for sh in shells:
        if sh not in blocks_by_shell:
            continue
        blk = dict(blocks_by_shell[sh])
        blk.setdefault("C_term_one_inf_hi", blk.get("C_term_hi"))
        blocks.append(blk)
    D = blocks[0]["D"] if blocks else 0
    full = blocks[0] if blocks else {}
    manifest = {
        "n": n,
        "D_full": D,
        "n_shells": len(shells),
        "radii_count": len(radii),
        "bound_method": "one_inf_only",
        "rational_coefficients": rational,
        "repair_tag": (
            "C0008_cert_repair_full_one_inf_rational"
            if rational
            else "C0008_cert_repair_full_one_inf"
        ),
        "repair_complete": complete,
        "blocks": blocks,
        "full_at_T_one_inf": full,
        "full_at_T_best": full,
        "arithmetic_backend": ArithmeticBackend(precision_bits=prec).as_dict(),
        "convention": (
            "M_r: output on shell r only. Repair baseline one_inf_only per-shell pass. "
            "Route A' cluster assembly uses L-0073R residual reweight."
            + (
                " Rational mpfr_leray + mpfr_vdot_real coefficients."
                if rational
                else ""
            )
        ),
    }
    payload = dict(manifest)
    payload.pop("manifest_sha256", None)
    manifest["manifest_sha256"] = hashlib.sha256(
        json.dumps(payload, sort_keys=True, default=str).encode()
    ).hexdigest()
    return manifest


def save_repair_manifest(manifest: dict, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return path


def repair_progress(n: int = 24, *, rational: bool = False) -> dict:
    radii = tuple(all_dealias_radii(n))
    ck_path, partial_path, hb_path, _ = _paths(n, rational=rational)
    done = _load_done(n=n, radii=radii, prec=128, resume=True, rational=rational)
    hb = hb_path.read_text(encoding="utf-8").strip() if hb_path.is_file() else ""
    return {
        "n": n,
        "rational_coefficients": rational,
        "n_total": len(radii),
        "n_completed": len(done),
        "complete": len(done) >= len(radii),
        "checkpoint": str(ck_path),
        "partial_manifest": str(partial_path) if partial_path.is_file() else None,
        "last_heartbeat": hb,
    }


def build_repair_full_one_inf_manifest(
    n: int = 24,
    *,
    prec: int = 128,
    workers: int = 1,
    resume: bool = True,
    shell_limit: int | None = None,
    out_path: Path | None = None,
    rational_coefficients: bool = False,
) -> dict:
    ck_path, partial_path, hb_path, ck_shell_dir = _paths(n, rational=rational_coefficients)
    radii = tuple(all_dealias_radii(n))
    done = _load_done(
        n=n, radii=radii, prec=prec, resume=resume, rational=rational_coefficients
    )
    pending = [sh for sh in radii if sh not in done]
    if shell_limit is not None:
        pending = pending[:shell_limit]

    n_total = len(radii)
    print(
        f"repair one_inf: n={n} shells={n_total} done={len(done)} pending={len(pending)} "
        f"rational={rational_coefficients}",
        flush=True,
    )
    _touch_heartbeat(hb_path, f"batch n={n} done={len(done)} pending={len(pending)}")

    if pending:
        ck_shell_dir.mkdir(parents=True, exist_ok=True)
        ck_dir = str(ck_shell_dir)

        def _finish(sh: int, res: dict) -> None:
            done[sh] = res
            _save_ck(
                n=n, radii=radii, prec=prec, done=done, rational=rational_coefficients
            )
            partial = assemble_repair_manifest(
                n=n,
                radii=radii,
                blocks_by_shell=done,
                prec=prec,
                complete=len(done) >= n_total,
                rational=rational_coefficients,
            )
            save_repair_manifest(partial, partial_path)
            _touch_heartbeat(hb_path, f"shell {sh} done [{len(done)}/{n_total}]")
            print(f"  shell {sh} C_term_hi={res['C_term_hi']} [{len(done)}/{n_total}]", flush=True)

        if workers > 1:
            while True:
                pending = [sh for sh in radii if sh not in done]
                if shell_limit is not None:
                    pending = pending[: max(0, shell_limit - len(done))]
                if not pending:
                    break
                effective = pick_workers(workers, n_pending=len(pending))
                batch = pending[:effective]
                print(
                    f"  batch: shells {batch} workers={effective} "
                    f"avail_ram={_available_ram_gb():.1f}GB [{len(done)}/{n_total}]",
                    flush=True,
                )
                tasks = [
                    (n, list(radii), sh, prec, ck_dir, str(hb_path), rational_coefficients)
                    for sh in batch
                ]
                failed: list[tuple[int, BaseException]] = []
                with ProcessPoolExecutor(max_workers=len(batch)) as ex:
                    futs = {ex.submit(_worker_one_shell, t): t[2] for t in tasks}
                    for fut in as_completed(futs):
                        sh = futs[fut]
                        try:
                            _finish(sh, fut.result())
                        except Exception as exc:
                            failed.append((sh, exc))
                            print(f"  shell {sh} worker failed: {exc!r}", flush=True)
                            traceback.print_exc()
                for sh, exc in failed:
                    print(f"  retry shell {sh} sequentially after worker failure...", flush=True)
                    try:
                        ck_path_sh = _resolve_shell_ck_path(
                            n, sh, ck_shell_dir, rational=rational_coefficients
                        )
                        bnd = certified_L_op_one_inf_single_shell_hi(
                            n,
                            radii,
                            sh,
                            prec=prec,
                            checkpoint_path=ck_path_sh,
                            heartbeat_path=hb_path,
                            rational_coefficients=rational_coefficients,
                        )
                        _finish(sh, {**bnd.as_dict(), "C_term_one_inf_hi": bnd.C_term_hi})
                    except Exception:
                        raise
        else:
            for sh in pending:
                print(f"  computing one_inf shell {sh}...", flush=True)
                ck_path_sh = _resolve_shell_ck_path(
                    n, sh, ck_shell_dir, rational=rational_coefficients
                )
                bnd = certified_L_op_one_inf_single_shell_hi(
                    n,
                    radii,
                    sh,
                    prec=prec,
                    checkpoint_path=ck_path_sh,
                    heartbeat_path=hb_path,
                    rational_coefficients=rational_coefficients,
                )
                _finish(sh, {**bnd.as_dict(), "C_term_one_inf_hi": bnd.C_term_hi})

    complete = len(done) >= n_total
    manifest = assemble_repair_manifest(
        n=n,
        radii=radii,
        blocks_by_shell=done,
        prec=prec,
        complete=complete,
        rational=rational_coefficients,
    )
    if complete:
        if out_path is None:
            tag = _rational_tag(rational_coefficients)
            out = data_root() / f"shell_manifest_repair_full_one_inf_n{n}{tag}.json"
        else:
            out = out_path
    else:
        out = out_path or partial_path
    save_repair_manifest(manifest, out)
    print(f"saved {out} complete={complete} sha={manifest['manifest_sha256'][:16]}...", flush=True)
    return manifest


def run_watchdog(
    n: int = 24,
    *,
    prec: int = 128,
    workers: int = 1,
    retry_delay: float = 10.0,
    max_retries: int = 5000,
    adaptive: bool = False,
    min_workers: int = 1,
    max_workers: int | None = None,
    rational_coefficients: bool = False,
) -> dict:
    """Auto-restart repair pass until all shells complete.

    With ``adaptive=True``, reduce parallel workers after OOM/crash and persist state.
    Intra-shell checkpoints survive worker/batch failures.
    """
    radii = tuple(all_dealias_radii(n))
    n_total = len(radii)
    state = _load_watchdog_state(n, rational=rational_coefficients)
    cur_workers = int(state.get("workers", workers))
    if max_workers is not None:
        cur_workers = min(cur_workers, max_workers)
    cur_workers = max(min_workers, cur_workers)
    last_exc: BaseException | None = None

    for attempt in range(1, max_retries + 1):
        prog = repair_progress(n, rational=rational_coefficients)
        if prog["n_completed"] >= n_total:
            print(f"watchdog: all {n_total} shells complete.", flush=True)
            _save_watchdog_state(
                n,
                {
                    "workers": cur_workers,
                    "status": "complete",
                    "n_completed": n_total,
                    "rational_coefficients": rational_coefficients,
                    "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                },
                rational=rational_coefficients,
            )
            return build_repair_full_one_inf_manifest(
                n,
                prec=prec,
                workers=min_workers,
                resume=True,
                rational_coefficients=rational_coefficients,
            )

        pending_n = n_total - prog["n_completed"]
        if adaptive:
            cur_workers = pick_workers(workers, n_pending=pending_n)
            cur_workers = max(min_workers, min(cur_workers, workers))
            if max_workers is not None:
                cur_workers = min(cur_workers, max_workers)
        else:
            cur_workers = workers

        print(
            f"watchdog: attempt {attempt} ({prog['n_completed']}/{n_total}) "
            f"workers={cur_workers} avail_ram={_available_ram_gb():.1f}GB...",
            flush=True,
        )
        _save_watchdog_state(
            n,
            {
                "workers": cur_workers,
                "status": "running",
                "attempt": attempt,
                "n_completed": prog["n_completed"],
                "rational_coefficients": rational_coefficients,
                "avail_ram_gb": round(_available_ram_gb(), 2),
                "last_error": repr(last_exc) if last_exc else None,
                "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            },
            rational=rational_coefficients,
        )
        last_exc = None
        try:
            build_repair_full_one_inf_manifest(
                n,
                prec=prec,
                workers=cur_workers,
                resume=True,
                shell_limit=None,
                rational_coefficients=rational_coefficients,
            )
        except KeyboardInterrupt:
            raise
        except Exception as exc:
            last_exc = exc
            print(f"watchdog: crash on attempt {attempt}:", flush=True)
            traceback.print_exc()
            if adaptive:
                new_w = _backoff_workers(cur_workers, exc, min_workers=min_workers)
                if new_w < cur_workers:
                    print(
                        f"watchdog: backing off workers {cur_workers} -> {new_w}",
                        flush=True,
                    )
                    cur_workers = new_w
        time.sleep(retry_delay)
    raise RuntimeError(f"watchdog: exceeded {max_retries} retries")
