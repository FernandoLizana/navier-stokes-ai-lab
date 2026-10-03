"""
External supervisor for C-0008 rational full-dealias regen.

Survives Cursor/IDE terminal death: run via start_c0008_rational_supervisor.ps1
(in external PowerShell), not as a Cursor background shell job.

- Restarts child regen on crash or exit before completion.
- Kills and restarts if heartbeat file is stale (stuck worker).
- Stop cleanly: scripts/stop_c0008_rational_supervisor.ps1
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
PYTHON = Path(sys.executable)


def _tag(n: int, rational: bool) -> str:
    return f"n{n}{'_rational' if rational else ''}"


def _paths(n: int, rational: bool) -> dict[str, Path]:
    t = _tag(n, rational)
    base = REPO / "experiments" / "terminal_weighted"
    return {
        "stop_flag": base / f"supervisor_stop_{t}.flag",
        "state": base / f"supervisor_state_{t}.json",
        "log": base / f"supervisor_{t}.log",
        "child_pid": base / f"supervisor_child_pid_{t}.txt",
        "heartbeat": base / f"repair_one_inf_heartbeat_{t}.txt",
    }


def _log(path: Path, msg: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    line = f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')} {msg}\n"
    with path.open("a", encoding="utf-8") as f:
        f.write(line)
    print(line, end="", flush=True)


def _save_state(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def _heartbeat_age_sec(hb: Path) -> float | None:
    if not hb.is_file():
        return None
    return max(0.0, time.time() - hb.stat().st_mtime)


def _progress(n: int, rational: bool) -> dict:
    sys.path.insert(0, str(REPO))
    from ns_exploration.terminal_weighted.repair_full_manifest import repair_progress

    return repair_progress(n, rational=rational)


def _kill_tree(pid: int) -> None:
    if sys.platform == "win32":
        subprocess.run(
            ["taskkill", "/F", "/T", "/PID", str(pid)],
            capture_output=True,
            check=False,
        )
    else:
        import os
        import signal

        try:
            os.killpg(os.getpgid(pid), signal.SIGTERM)
        except ProcessLookupError:
            pass


def _child_cmd(args: argparse.Namespace) -> list[str]:
    cmd = [
        str(PYTHON),
        str(REPO / "scripts" / "regenerate_c0008_full_dealias_one_inf.py"),
        "--n",
        str(args.n),
        "--workers",
        str(args.workers),
        "--watchdog",
        "--retry-delay",
        str(args.retry_delay),
    ]
    if args.rational:
        cmd.append("--rational")
    if args.adaptive:
        cmd.append("--adaptive")
    if args.max_workers is not None:
        cmd.extend(["--max-workers", str(args.max_workers)])
    if args.min_workers != 1:
        cmd.extend(["--min-workers", str(args.min_workers)])
    return cmd


def supervise(args: argparse.Namespace) -> int:
    paths = _paths(args.n, args.rational)
    if paths["stop_flag"].is_file():
        paths["stop_flag"].unlink()

    attempt = 0
    while not paths["stop_flag"].is_file():
        prog = _progress(args.n, args.rational)
        if prog.get("complete"):
            _log(paths["log"], f"COMPLETE n_completed={prog['n_completed']}/{prog['n_total']}")
            _save_state(
                paths["state"],
                {
                    "status": "complete",
                    "n_completed": prog["n_completed"],
                    "n_total": prog["n_total"],
                    "updated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                },
            )
            return 0

        attempt += 1
        cmd = _child_cmd(args)
        _log(
            paths["log"],
            f"attempt={attempt} spawn {' '.join(cmd)} "
            f"progress={prog['n_completed']}/{prog['n_total']}",
        )
        child_log = paths["log"].with_name(paths["log"].stem + "_child.log")
        with child_log.open("a", encoding="utf-8") as lf:
            lf.write(
                f"\n--- attempt {attempt} {datetime.now(timezone.utc).isoformat()} ---\n"
            )
            lf.flush()
            proc = subprocess.Popen(
                cmd,
                cwd=str(REPO),
                stdout=lf,
                stderr=subprocess.STDOUT,
                text=True,
            )
            paths["child_pid"].write_text(str(proc.pid), encoding="utf-8")
            _save_state(
                paths["state"],
                {
                    "status": "running",
                    "supervisor_pid": __import__("os").getpid(),
                    "child_pid": proc.pid,
                    "attempt": attempt,
                    "n_completed": prog["n_completed"],
                    "n_total": prog["n_total"],
                    "stale_minutes": args.stale_minutes,
                    "updated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
                },
            )

            stale_killed = False
            while proc.poll() is None:
                if paths["stop_flag"].is_file():
                    _log(paths["log"], f"stop flag seen; terminating child pid={proc.pid}")
                    _kill_tree(proc.pid)
                    proc.wait(timeout=120)
                    paths["child_pid"].unlink(missing_ok=True)
                    _save_state(
                        paths["state"],
                        {
                            "status": "stopped",
                            "updated_at": datetime.now(timezone.utc).strftime(
                                "%Y-%m-%dT%H:%M:%SZ"
                            ),
                        },
                    )
                    return 0

                age = _heartbeat_age_sec(paths["heartbeat"])
                if age is not None and age > args.stale_minutes * 60:
                    _log(
                        paths["log"],
                        f"STALE heartbeat {age:.0f}s > {args.stale_minutes * 60}s; "
                        f"killing child pid={proc.pid}",
                    )
                    _kill_tree(proc.pid)
                    proc.wait(timeout=120)
                    stale_killed = True
                    break

                time.sleep(args.poll_seconds)

        rc = proc.returncode if proc.returncode is not None else proc.wait()
        paths["child_pid"].unlink(missing_ok=True)

        prog = _progress(args.n, args.rational)
        if prog.get("complete"):
            _log(paths["log"], f"COMPLETE after child exit rc={rc}")
            _save_state(
                paths["state"],
                {"status": "complete", "n_completed": prog["n_completed"], "n_total": prog["n_total"]},
            )
            return 0

        reason = "stale_heartbeat" if stale_killed else f"exit rc={rc}"
        _log(
            paths["log"],
            f"child ended ({reason}) progress={prog['n_completed']}/{prog['n_total']}; "
            f"restart in {args.restart_delay}s",
        )
        _save_state(
            paths["state"],
            {
                "status": "restarting",
                "last_reason": reason,
                "n_completed": prog["n_completed"],
                "attempt": attempt,
                "updated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            },
        )
        time.sleep(args.restart_delay)

    return 0


def main() -> int:
    p = argparse.ArgumentParser(description="Supervise C-0008 rational regen (auto-restart)")
    p.add_argument("--n", type=int, default=24)
    p.add_argument("--rational", action="store_true", default=True)
    p.add_argument("--workers", type=int, default=8)
    p.add_argument("--max-workers", type=int, default=None)
    p.add_argument("--min-workers", type=int, default=1)
    p.add_argument("--adaptive", action="store_true", default=True)
    p.add_argument("--stale-minutes", type=float, default=25.0, help="kill/restart if heartbeat older")
    p.add_argument("--poll-seconds", type=float, default=30.0)
    p.add_argument("--restart-delay", type=float, default=15.0)
    p.add_argument("--retry-delay", type=float, default=10.0, help="passed to inner watchdog")
    args = p.parse_args()
    return supervise(args)


if __name__ == "__main__":
    raise SystemExit(main())
