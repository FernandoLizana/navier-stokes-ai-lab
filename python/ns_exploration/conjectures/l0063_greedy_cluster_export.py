"""
L-0063: Greedy C-R-0012 sparse export status (cluster prep).

Runs or records export of greedy 41-shell support for cluster SOS.
Does not launch cluster SDP on laptop.

FINITE Galerkin only. Not Clay.
"""

from __future__ import annotations

import json
import subprocess
import sys
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0047_greedy_fullsym import SUPPORT_CR0012
from ns_exploration.conjectures.l0055_greedy_sos_feasibility import lemma_l0055
from ns_exploration.conjectures.l0057_greedy_cluster_sos import (
    GREEDY_DATA,
    lemma_l0057,
    save_lemma_l0057,
    write_cluster_runbook,
    write_slurm_template,
)

REPO = Path(__file__).resolve().parents[3]


@dataclass
class GalerkinBoundL0063:
    lemma_id: str = "L-0063"
    route: str = "B"
    status: str = "cluster_export_pending"
    evidence_level: str = "N7"
    n: int = 24
    n_shells: int = 41
    greedy_D: int = 1772
    export_dir: str = ""
    export_complete: bool = False
    nnz_M: int = 0
    C_ub_sos_extrap: float = 0.0
    cluster_runbook: str = ""
    clay_implication: str = "None. Greedy export prep only; not Clay."
    notes: str = ""

    def as_dict(self) -> dict:
        return asdict(self)


def _load_meta() -> dict:
    p = GREEDY_DATA / "meta.json"
    if not p.is_file():
        return {}
    return json.loads(p.read_text(encoding="utf-8"))


def run_greedy_export(n: int = 24, timeout_s: int = 7200) -> dict:
    """Export sparse fullsym M for greedy C-R-0012 (may take tens of minutes)."""
    GREEDY_DATA.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        "-m",
        "ns_exploration.experiments.export_sos_cubic",
        "--n",
        str(n),
        "--cr0012",
        "--out",
        str(GREEDY_DATA),
    ]
    env = {"PYTHONPATH": str(REPO / "python")}
    r = subprocess.run(
        cmd,
        cwd=str(REPO),
        capture_output=True,
        text=True,
        timeout=timeout_s,
        env={**dict(**__import__("os").environ), **env},
    )
    return {"returncode": r.returncode, "stdout": r.stdout[-2000:], "stderr": r.stderr[-2000:]}


def lemma_l0063(n: int = 24) -> GalerkinBoundL0063:
    l55 = lemma_l0055(n=n)
    meta = _load_meta()
    complete = bool(meta.get("D")) and (GREEDY_DATA / "M_csr.npz").is_file()
    status = "cluster_export_ready" if complete else "cluster_export_pending"
    notes = (
        f"Greedy C-R-0012 export: {len(SUPPORT_CR0012)} shells, D≈{l55.greedy_D}. "
        f"Extrap SOS C_ub≈{l55.C_ub_sos_extrap:.3f}. "
        f"Export {'done' if complete else 'pending'} → cluster SLURM (L-0057)."
    )
    return GalerkinBoundL0063(
        n=n,
        greedy_D=l55.greedy_D,
        export_dir=str(GREEDY_DATA),
        export_complete=complete,
        nnz_M=int(meta.get("nnz_M", 0)),
        C_ub_sos_extrap=l55.C_ub_sos_extrap,
        cluster_runbook="tools/sos_julia/cluster/greedy_cr0012/README.md",
        status=status,
        notes=notes,
    )


def save_lemma_l0063(
    bound: GalerkinBoundL0063,
    path: str | Path = "conjectures/active/L-0063.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
