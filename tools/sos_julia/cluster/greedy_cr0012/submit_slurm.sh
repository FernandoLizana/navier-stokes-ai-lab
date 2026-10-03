#!/bin/bash
#SBATCH --job-name=greedy-cr0012-sos
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=120:00:00
#SBATCH --output=greedy_cr0012_%j.log

set -euo pipefail
REPO="${REPO:-$(pwd)}"
cd "$REPO"
export PYTHONPATH=python

# Stage 1: export (skip if meta.json exists)
if [[ ! -f tools/sos_julia/data/greedy_cr0012/meta.json ]]; then
  python -m ns_exploration.experiments.export_sos_cubic --n 24 --cr0012 \
    --out tools/sos_julia/data/greedy_cr0012
fi

# Stage 2: TSSOS COSMO
cd tools/sos_julia
julia --project=tssos_env run_tssos_cosmo.jl data/greedy_cr0012

# Stage 3: Gram export (optional; large disk)
julia --project=tssos_env export_tssos_gram.jl data/greedy_cr0012

cd "$REPO"
python -m ns_exploration.experiments.sprint02_l0057_cluster
