"""
L-0057: Cluster job pack for greedy C-R-0012 SOS (D≈1772).

Packages export metadata, resource estimates (L-0055), and cluster runbook.
Does NOT launch the cluster SDP locally.

FINITE Galerkin only. Not Clay.
"""

from __future__ import annotations

import json
from dataclasses import asdict, dataclass
from pathlib import Path

from ns_exploration.conjectures.l0027_low_slab_cubic import lemma_l0027
from ns_exploration.conjectures.l0047_greedy_fullsym import SUPPORT_CR0012
from ns_exploration.conjectures.l0051_structured_ladder import C_GREEDY_41_SHELLS
from ns_exploration.conjectures.l0055_greedy_sos_feasibility import lemma_l0055

GREEDY_DATA = Path("tools/sos_julia/data/greedy_cr0012")
CLUSTER_DIR = Path("tools/sos_julia/cluster/greedy_cr0012")


@dataclass
class ClusterStage:
    id: str
    command: str
    est_hours: float
    est_ram_gb: float
    notes: str


@dataclass
class GalerkinBoundL0057:
    lemma_id: str = "L-0057"
    route: str = "B"
    status: str = "cluster_ready"
    evidence_level: str = "N7"
    n: int = 24
    C_dagger: float = 0.0
    C_shor: float = C_GREEDY_41_SHELLS
    greedy_D: int = 0
    n_shells: int = 41
    C_ub_sos_extrap: float = 0.0
    est_sdp_variables_M: float = 0.0
    est_solve_hours: float = 0.0
    est_gram_export_gb: float = 0.0
    export_dir: str = ""
    cluster_dir: str = ""
    stages: list[ClusterStage] | None = None
    clay_implication: str = (
        "None. Cluster job pack only; not a SOS certificate; not Clay."
    )
    notes: str = ""

    def as_dict(self) -> dict:
        d = asdict(self)
        if self.stages:
            d["stages"] = [asdict(s) for s in self.stages]
        return d


def _load_greedy_meta() -> dict:
    meta_path = GREEDY_DATA / "meta.json"
    if not meta_path.is_file():
        return {}
    return json.loads(meta_path.read_text(encoding="utf-8"))


def lemma_l0057(n: int = 24) -> GalerkinBoundL0057:
    Cd = lemma_l0027(n=n, empirical=False).C_dagger
    l55 = lemma_l0055(n=n)
    meta = _load_greedy_meta()
    D = int(meta.get("D", l55.greedy_D))
    stages = [
        ClusterStage(
            id="export",
            command=(
                "PYTHONPATH=python python -m ns_exploration.experiments.export_sos_cubic "
                f"--n {n} --cr0012 --out {GREEDY_DATA.as_posix()}"
            ),
            est_hours=0.5,
            est_ram_gb=16.0,
            notes="Sparse fullsym M export for 41-shell greedy witness.",
        ),
        ClusterStage(
            id="tssos_cosmo",
            command=(
                f"cd tools/sos_julia && julia --project=tssos_env "
                f"run_tssos_cosmo.jl {GREEDY_DATA.as_posix()}"
            ),
            est_hours=l55.est_solve_hours,
            est_ram_gb=64.0,
            notes=f"Order-2 TSSOS TS=MD COSMO; est ~{l55.est_sdp_variables/1e6:.1f}M SDP vars.",
        ),
        ClusterStage(
            id="gram_export",
            command=(
                f"cd tools/sos_julia && julia --project=tssos_env "
                f"export_tssos_gram.jl {GREEDY_DATA.as_posix()}"
            ),
            est_hours=8.0,
            est_ram_gb=128.0,
            notes=f"Rational Gram export; est ~{l55.est_gram_export_gb:.0f} GB JSON.",
        ),
        ClusterStage(
            id="n5_verify",
            command="PYTHONPATH=python python -m ns_exploration.experiments.sprint02_l0057_cluster",
            est_hours=2.0,
            est_ram_gb=32.0,
            notes="Python N5 Gram verify + certificate (post-cluster).",
        ),
    ]
    notes = (
        f"Greedy C-R-0012 cluster pack: D={D}, C_Shor≈{C_GREEDY_41_SHELLS:.4f}. "
        f"Extrap SOS C_ub≈{l55.C_ub_sos_extrap:.3f} (< C_dagger≈{Cd:.3f}). "
        f"Est SDP ~{l55.est_sdp_variables/1e6:.1f}M vars, ~{l55.est_solve_hours:.0f}h solve, "
        f"Gram ~{l55.est_gram_export_gb:.0f}GB. Laptop infeasible; submit via SLURM template."
    )
    return GalerkinBoundL0057(
        n=n,
        C_dagger=Cd,
        greedy_D=D,
        n_shells=len(SUPPORT_CR0012),
        C_ub_sos_extrap=l55.C_ub_sos_extrap,
        est_sdp_variables_M=l55.est_sdp_variables / 1e6,
        est_solve_hours=l55.est_solve_hours,
        est_gram_export_gb=l55.est_gram_export_gb,
        export_dir=str(GREEDY_DATA),
        cluster_dir=str(CLUSTER_DIR),
        stages=stages,
        notes=notes,
    )


def write_cluster_runbook(bound: GalerkinBoundL0057) -> Path:
    CLUSTER_DIR.mkdir(parents=True, exist_ok=True)
    md = [
        "# Greedy C-R-0012 SOS — cluster runbook (L-0057)",
        "",
        "> Finite Galerkin T³, N≤24. **Not continuum. Not Clay.**",
        "",
        f"- **Support:** {len(SUPPORT_CR0012)} shells, D≈{bound.greedy_D}",
        f"- **C_Shor:** {bound.C_shor:.6f}",
        f"- **Extrap SOS C_ub:** {bound.C_ub_sos_extrap:.4f} (< C_dagger {bound.C_dagger:.6f})",
        f"- **Est SDP vars:** ~{bound.est_sdp_variables_M:.1f}M",
        f"- **Est solve:** ~{bound.est_solve_hours:.0f} h",
        f"- **Est Gram export:** ~{bound.est_gram_export_gb:.0f} GB",
        "",
        "## Stages",
        "",
    ]
    for s in bound.stages or []:
        md.extend(
            [
                f"### {s.id}",
                f"- RAM: ~{s.est_ram_gb:.0f} GB",
                f"- Time: ~{s.est_hours:.1f} h",
                f"- Notes: {s.notes}",
                "",
                "```bash",
                s.command,
                "```",
                "",
            ]
        )
    path = CLUSTER_DIR / "README.md"
    path.write_text("\n".join(md), encoding="utf-8")
    return path


def write_slurm_template(bound: GalerkinBoundL0057) -> Path:
    CLUSTER_DIR.mkdir(parents=True, exist_ok=True)
    script = f"""#!/bin/bash
#SBATCH --job-name=greedy-cr0012-sos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --output=greedy_cr0012_%j.log

set -euo pipefail
REPO="${{REPO:-$(pwd)}}"
cd "$REPO"
export PYTHONPATH=python

# Stage 1: export (skip if meta.json exists)
if [[ ! -f tools/sos_julia/data/greedy_cr0012/meta.json ]]; then
  python -m ns_exploration.experiments.export_sos_cubic --n {bound.n} --cr0012 \\
    --out tools/sos_julia/data/greedy_cr0012
fi

# Stage 2: TSSOS COSMO
cd tools/sos_julia
julia --project=tssos_env run_tssos_cosmo.jl data/greedy_cr0012

# Stage 3: Gram export (optional; large disk)
julia --project=tssos_env export_tssos_gram.jl data/greedy_cr0012

cd "$REPO"
python -m ns_exploration.experiments.sprint02_l0057_cluster
"""
    path = CLUSTER_DIR / "submit_slurm.sh"
    path.write_text(script, encoding="utf-8")
    return path


def save_lemma_l0057(
    bound: GalerkinBoundL0057,
    path: str | Path = "conjectures/active/L-0057.json",
) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(bound.as_dict(), indent=2), encoding="utf-8")
