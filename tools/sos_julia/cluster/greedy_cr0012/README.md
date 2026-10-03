# Greedy C-R-0012 SOS — cluster runbook (L-0057)

> Finite Galerkin T³, N≤24. **Not continuum. Not Clay.**

- **Support:** 41 shells, D≈1772
- **C_Shor:** 9.503700
- **Extrap SOS C_ub:** 5.3886 (< C_dagger 9.562010)
- **Est SDP vars:** ~8.6M
- **Est solve:** ~113 h
- **Est Gram export:** ~4353 GB

## Stages

### export
- RAM: ~16 GB
- Time: ~0.5 h
- Notes: Sparse fullsym M export for 41-shell greedy witness.

```bash
PYTHONPATH=python python -m ns_exploration.experiments.export_sos_cubic --n 24 --cr0012 --out tools/sos_julia/data/greedy_cr0012
```

### tssos_cosmo
- RAM: ~64 GB
- Time: ~113.2 h
- Notes: Order-2 TSSOS TS=MD COSMO; est ~8.6M SDP vars.

```bash
cd tools/sos_julia && julia --project=tssos_env run_tssos_cosmo.jl tools/sos_julia/data/greedy_cr0012
```

### gram_export
- RAM: ~128 GB
- Time: ~8.0 h
- Notes: Rational Gram export; est ~4353 GB JSON.

```bash
cd tools/sos_julia && julia --project=tssos_env export_tssos_gram.jl tools/sos_julia/data/greedy_cr0012
```

### n5_verify
- RAM: ~32 GB
- Time: ~2.0 h
- Notes: Python N5 Gram verify + certificate (post-cluster).

```bash
PYTHONPATH=python python -m ns_exploration.experiments.sprint02_l0057_cluster
```
