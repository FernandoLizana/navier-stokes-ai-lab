#!/usr/bin/env python3
"""
Portable job supervisor for C-0008 CPU compute (Ubuntu / Linux / macOS / Windows).

Reads compute_pack/jobs.json. All checkpoints/logs under compute_pack/data/.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

PACK_ROOT = Path(__file__).resolve().parent
REPO_ROOT = PACK_ROOT.parent
PYTHON_DIR = REPO_ROOT / "python"
JOBS_FILE = PACK_ROOT / "jobs.json"


def _load_jobs() -> list[dict]:
    data = json.loads(JOBS_FILE.read_text(encoding="utf-8"))
    return data["jobs"]


def _job_by_id(job_id: str) -> dict:
    for j in _load_jobs():
        if j["id"] == job_id:
            return j
    raise SystemExit(f"unknown job_id: {job_id}")


def _paths(job_id: str) -> dict[str, Path]:
    data = PACK_ROOT / "data" / "terminal_weighted"
    logs = PACK_ROOT / "data" / "logs"
    data.mkdir(parents=True, exist_ok=True)
    logs.mkdir(parents=True, exist_ok=True)
    return {
        "stop_flag": logs / f"stop_{job_id}.flag",
        "supervisor_log": logs / f"supervisor_{job_id}.log",
        "child_log": logs / f"child_{job_id}.log",
        "state": logs / f"state_{job_id}.json",
        "child_pid": logs / f"child_pid_{job_id}.txt",
        "data": data,
    }


def _setup_env() -> dict[str, str]:
    env = os.environ.copy()
    env["C0008_DATA_DIR"] = str(PACK_ROOT / "data" / "terminal_weighted")
    env["C0008_LOG_DIR"] = str(PACK_ROOT / "data" / "logs")
    env["C0008_CERT_DIR"] = str(PACK_ROOT / "data" / "certificates")
    env["PYTHONPATH"] = str(PYTHON_DIR)
    return env


def _log(path: Path, msg: str) -> None:
    line = f"{datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')} {msg}\n"
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("a", encoding="utf-8") as f:
        f.write(line)
    print(line, end="", flush=True)


def _heartbeat_age(data_dir: Path, hb_name: str | None) -> float | None:
    if not hb_name:
        return None
    hb = data_dir / hb_name
    if not hb.is_file():
        return None
    return max(0.0, time.time() - hb.stat().st_mtime)


def _is_complete(job: dict) -> bool:
    chk = job.get("complete_check")
    if not chk:
        return False
    if chk.get("type") != "repair_progress":
        return False
    sys.path.insert(0, str(PYTHON_DIR))
    os.environ.setdefault("C0008_DATA_DIR", str(PACK_ROOT / "data" / "terminal_weighted"))
    from ns_exploration.terminal_weighted.repair_full_manifest import repair_progress

    prog = repair_progress(chk["n"], rational=chk.get("rational", False))
    return bool(prog.get("complete"))


def _deps_complete(job: dict, done: set[str]) -> bool:
    return all(d in done for d in job.get("depends_on", []))


def _kill_tree(pid: int) -> None:
    if sys.platform == "win32":
        subprocess.run(["taskkill", "/F", "/T", "/PID", str(pid)], check=False, capture_output=True)
    else:
        import signal

        try:
            os.killpg(os.getpgid(pid), signal.SIGTERM)
        except (ProcessLookupError, PermissionError):
            try:
                os.kill(pid, signal.SIGTERM)
            except ProcessLookupError:
                pass


def _resolve_workers(requested: int) -> int:
    if requested > 0:
        return requested
    sys.path.insert(0, str(PACK_ROOT))
    from ram_workers import recommend_workers

    return recommend_workers(requested_max=int(os.environ.get("C0008_WORKERS_MAX", "32")))


def _run_cmd(cmd: list[str], env: dict, log_path: Path) -> int:
    with log_path.open("a", encoding="utf-8") as lf:
        lf.write(f"\n--- {datetime.now(timezone.utc).isoformat()} {' '.join(cmd)} ---\n")
        lf.flush()
        proc = subprocess.Popen(
            [sys.executable, *cmd],
            cwd=str(PYTHON_DIR),
            env=env,
            stdout=lf,
            stderr=subprocess.STDOUT,
            text=True,
            start_new_session=sys.platform != "win32",
        )
        return proc.wait()


def supervise_job(job_id: str, *, workers: int, adaptive: bool) -> int:
    workers = _resolve_workers(workers)
    job = _job_by_id(job_id)
    paths = _paths(job_id)
    env = _setup_env()
    if paths["stop_flag"].is_file():
        paths["stop_flag"].unlink()

    stale_min = float(job.get("stale_minutes", 25))
    poll = 30.0
    restart_delay = 15.0
    hb_name = job.get("heartbeat")
    attempt = 0

    base_cmd = list(job["command"])
    if "--workers" not in base_cmd and job.get("workers_default"):
        base_cmd.extend(["--workers", str(workers)])
    if adaptive and "--adaptive" not in base_cmd:
        base_cmd.append("--adaptive")

    while not paths["stop_flag"].is_file():
        if _is_complete(job):
            _log(paths["supervisor_log"], f"COMPLETE job={job_id}")
            break

        attempt += 1
        _log(paths["supervisor_log"], f"attempt={attempt} spawn {base_cmd}")
        with paths["child_log"].open("a", encoding="utf-8") as lf:
            lf.write(f"\n=== attempt {attempt} ===\n")
            proc = subprocess.Popen(
                [sys.executable, *base_cmd],
                cwd=str(PYTHON_DIR),
                env=env,
                stdout=lf,
                stderr=subprocess.STDOUT,
                text=True,
                start_new_session=sys.platform != "win32",
            )
        paths["child_pid"].write_text(str(proc.pid), encoding="utf-8")

        stale_killed = False
        while proc.poll() is None:
            if paths["stop_flag"].is_file():
                _log(paths["supervisor_log"], "stop flag; terminating")
                _kill_tree(proc.pid)
                proc.wait(timeout=120)
                paths["child_pid"].unlink(missing_ok=True)
                return 0
            age = _heartbeat_age(paths["data"], hb_name)
            if age is not None and age > stale_min * 60:
                _log(paths["supervisor_log"], f"STALE heartbeat {age:.0f}s; kill pid={proc.pid}")
                _kill_tree(proc.pid)
                proc.wait(timeout=120)
                stale_killed = True
                break
            time.sleep(poll)

        rc = proc.returncode if proc.returncode is not None else proc.wait()
        paths["child_pid"].unlink(missing_ok=True)

        if _is_complete(job):
            _log(paths["supervisor_log"], f"COMPLETE after rc={rc}")
            break

        reason = "stale" if stale_killed else f"rc={rc}"
        _log(paths["supervisor_log"], f"child ended ({reason}); restart in {restart_delay}s")
        time.sleep(restart_delay)

    for post in job.get("post_commands", []):
        _log(paths["supervisor_log"], f"post: {post}")
        _run_cmd(list(post), env, paths["child_log"])

    paths["state"].write_text(
        json.dumps(
            {
                "job_id": job_id,
                "status": "complete" if _is_complete(job) else "stopped",
                "updated_at": datetime.now(timezone.utc).isoformat(),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    return 0


def run_queue(*, workers: int, adaptive: bool, only: str | None) -> int:
    workers = _resolve_workers(workers)
    jobs = sorted(_load_jobs(), key=lambda j: j.get("priority", 99))
    done: set[str] = set()
    for job in jobs:
        if only and job["id"] != only:
            continue
        if job.get("optional") and only is None:
            continue
        if not _deps_complete(job, done):
            _log(_paths(job["id"])["supervisor_log"], f"skip {job['id']}: deps pending")
            continue
        if _is_complete(job):
            done.add(job["id"])
            continue
        w = workers or job.get("workers_default", 4)
        supervise_job(job["id"], workers=w, adaptive=adaptive)
        if _is_complete(job):
            done.add(job["id"])
    return 0


def cmd_status(all_jobs: bool, job_id: str | None) -> None:
    env = _setup_env()
    os.environ.update({k: env[k] for k in ("C0008_DATA_DIR", "C0008_LOG_DIR")})
    sys.path.insert(0, str(PYTHON_DIR))
    from ns_exploration.terminal_weighted.repair_full_manifest import repair_progress

    ids = [job_id] if job_id else [j["id"] for j in _load_jobs()]
    for jid in ids:
        job = _job_by_id(jid)
        chk = job.get("complete_check")
        prog = None
        if chk and chk.get("type") == "repair_progress":
            prog = repair_progress(chk["n"], rational=chk.get("rational", False))
        print(json.dumps({"job_id": jid, "title": job["title"], "progress": prog, "complete": _is_complete(job)}, indent=2))


def cmd_stop(job_id: str | None) -> None:
    ids = [job_id] if job_id else [j["id"] for j in _load_jobs()]
    for jid in ids:
        p = _paths(jid)
        p["stop_flag"].touch()
        if p["child_pid"].is_file():
            try:
                _kill_tree(int(p["child_pid"].read_text().strip()))
            except ValueError:
                pass
    print("stop flags set")


def main() -> int:
    p = argparse.ArgumentParser(description="C-0008 portable compute supervisor")
    p.add_argument("job_id", nargs="?", help="job id from jobs.json (omit with --queue)")
    p.add_argument("--queue", action="store_true", help="run priority queue (non-optional jobs)")
    p.add_argument("--workers", type=int, default=0, help="0 = use ram_workers up to 32")
    p.add_argument("--adaptive", action="store_true", default=True)
    p.add_argument("--status", action="store_true")
    p.add_argument("--status-all", action="store_true")
    p.add_argument("--stop", action="store_true")
    p.add_argument("--list", action="store_true")
    args = p.parse_args()

    if args.list:
        for j in _load_jobs():
            opt = " [optional]" if j.get("optional") else ""
            print(f"{j['priority']:2d}  {j['id']:<28} ~{j['estimated_hours']}h  {j['title']}{opt}")
        return 0
    if args.status or args.status_all:
        cmd_status(args.status_all, args.job_id)
        return 0
    if args.stop:
        cmd_stop(args.job_id)
        return 0
    if args.queue:
        return run_queue(workers=args.workers, adaptive=args.adaptive, only=None)
    if not args.job_id:
        p.error("job_id required unless --queue/--list/--status")
    return supervise_job(args.job_id, workers=args.workers, adaptive=args.adaptive)


if __name__ == "__main__":
    raise SystemExit(main())
