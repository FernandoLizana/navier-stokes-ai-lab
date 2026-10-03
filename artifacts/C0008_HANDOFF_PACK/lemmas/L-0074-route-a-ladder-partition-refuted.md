# L-0074 — Route A ladder: partition family refuted

**Status:** **negative result** — finer cluster partitions do not improve L-0073

**Parent:** L-0073

**Evidence:** N5 full-dealias ladder (4 configs, recomputed from scratch)

## Result (2026-08-03)

| Config | clusters | Route A I_hi | best | Omega(T)^hi |
|--------|----------|--------------|------|-------------|
| **equal_12** | 11 | **37.92** | A | **54.23** |
| fine_low | 19 | 40.13 | A | 55.01 |
| equal_24 | 22 | 63.85 | A | 63.40 |
| equal_48 | 44 | 70.50 | D | 65.75 |

**Best remains L-0073 (`equal_12`).** Monotonic worsening with more equal clusters; `fine_low` second but still +0.8 Omega.

## Refutation scope

- Refutes: hope that **adaptive/finer cluster partitions** + min(Route A', Route D) beat 12 equal clusters on full dealias.
- Does **not** refute C-0008 or L-0073 certificate integrity.

## Artifacts

- Summary: `experiments/terminal_weighted/route_a_ladder_full.json`
- Report: `reports/C0008_ROUTE_A_LADDER_FULL.md`
- Script: `python/ns_exploration/experiments/sprint_c0008_route_a_ladder.py`

## Next viable paths

1. L-0027 stretch (all-IC C <= C_dagger)
2. Weighted triad / SOS subclasses
3. Non-partition terminal structure (new lemma class)
