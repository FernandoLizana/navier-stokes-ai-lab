"""Repository snapshot hashes for terminal-weighted research (spec §3)."""

from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_SNAPSHOT_FILES = [
    "docs/ROADMAP.md",
    "reports/PAQUETE_CHATGPT_REVIEW_COMPLETO.md",
    "reports/C0007_STRUCTURED_CLOSURE_MAP.md",
    "conjectures/active/C-0008.json",
    "conjectures/proved_restricted/C-0007.json",
    "conjectures/proved_restricted/L-0048.json",
    "certificates/CERT-L0070-all-ic-C0007-N24.json",
    "experiments/exploratory/sprint02_l0070/summary.json",
    "tools/sos_julia/data/scale_results.json",
    "python/ns_exploration/conjectures/l0048_fullsym_techo.py",
    "python/ns_exploration/conjectures/l0048_matrixfree_fullsym.py",
    "python/ns_exploration/conjectures/l0045_streaming_fullsym.py",
    "python/ns_exploration/experiments/sprintB_physical_tensor.py",
]


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def build_repo_snapshot(
    root: str | Path = ".",
    files: list[str] | None = None,
) -> dict:
    root = Path(root)
    files = files or DEFAULT_SNAPSHOT_FILES
    entries = []
    missing = []
    for rel in files:
        p = root / rel
        if not p.is_file():
            missing.append(rel)
            continue
        entries.append(
            {
                "path": rel.replace("\\", "/"),
                "sha256": sha256_file(p),
                "bytes": p.stat().st_size,
            }
        )
    return {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "root": str(root.resolve()),
        "files": entries,
        "missing": missing,
        "n_ok": len(entries),
    }


def write_repo_snapshot(
    path: str | Path = "experiments/terminal_weighted/repo_snapshot.json",
    root: str | Path = ".",
) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    snap = build_repo_snapshot(root=root)
    path.write_text(json.dumps(snap, indent=2), encoding="utf-8")
    return path
