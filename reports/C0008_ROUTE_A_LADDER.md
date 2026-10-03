# C-0008 Route A Ladder — Adaptive Clusters + min(Route A', Route D)

**Status:** full sweep **complete** (2026-08-03). Best: `equal_12` = L-0073 (I_hi=37.92, Omega=54.23). Verifier FAIL expected.

## Band pilot (shells 1–6, best bounds)

| Label | clusters | Route A | Route D | pick | Ω(T)^hi |
|-------|----------|---------|---------|------|---------|
| **equal_2** | 2 | **0.837** | 1.144 | A | **41.12** |
| equal_4 | 3 | 0.982 | 1.145 | A | 41.17 |
| equal_6 | 6 | 1.244 | 1.146 | D | 41.23 |
| fine_low | 6 | 1.244 | 1.146 | D | 41.23 |
| singleton | 6 | 1.244 | 1.146 | D | 41.23 |

**Best band:** `equal_2` Route A — Ω ≈ 41.12 (below M ≈ 41.284 on band only).

On band, finer partitions (singleton / fine_low) **lose** the cluster Frobenius gain; 2 coarse clusters remain optimal.

**Status:** full sweep **complete** (2026-08-03).

## Full-dealias results

| Label | clusters | Route A | Route D | pick | Omega |
|-------|----------|---------|---------|------|-------|
| **equal_12** | 11 | **37.92** | 62.48 | A | **54.23** |
| fine_low | 19 | 40.13 | 62.27 | A | 55.01 |
| equal_24 | 22 | 63.85 | 65.58 | A | 63.40 |
| equal_48 | 44 | 103.87 | 70.50 | D | 65.75 |

12 equal clusters optimal; finer partitions worsen bounds. L-0073 unchanged.

## Artifacts

- Full summary: `experiments/terminal_weighted/route_a_ladder_full.json`
- Full manifest: `experiments/terminal_weighted/shell_manifest_route_a_best_full.json`
- Finalize run: `experiments/terminal_weighted/route_a_ladder_finalize.json`
- Script: `python/ns_exploration/experiments/sprint_c0008_route_a_ladder.py`

## Method

- **Route A'**: per cluster `min(Σ B_r f(r), B_C max f(r))` with cluster+best bounds.
- **Route D (cluster-aware)**: `C(T)·T + Σ_C min(Σ B_r int_φ(r), B_C max int_φ)`.
- **Ladder**: adaptive partitions (`equal`, `fine_low`, `singleton`) × cluster counts.
