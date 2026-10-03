# L-0075 — Hybrid low-slab bound (triad + defect + SOS cubic)

**Status:** **best low-slab hybrid sketch** — does not close C-0008 sharp M

**Parent:** L-0027 / L-0030 / L-0031

**Evidence:** N7 ODE probes + subclass SOS certificates

## Method

On low slab Ω₀ ∈ [E₀, Ω★], combine:

- L-0030: min(hybrid Young/triad, spectral defect) enstrophy ODE
- L-0031: weighted triad + defect
- SOS subclass: cubic ODE with certified C_ub_hi (band123, greedy one-pol)

Take **minimum** worst-case Ω(T) over routes.

## Artifacts

- Summary: `experiments/terminal_weighted/hybrid_low_slab_summary.json`
- Report: `reports/C0008_HYBRID_LOW_SLAB.md`
- Script: `python/ns_exploration/experiments/sprint_c0008_hybrid_low_slab.py`

Run:

```powershell
python -m ns_exploration.experiments.sprint_c0008_hybrid_low_slab
```

**Note:** SOS routes are **subclass-only**. All-IC C ≤ C_† remains open (C_fullsym ≈ 25.9).
