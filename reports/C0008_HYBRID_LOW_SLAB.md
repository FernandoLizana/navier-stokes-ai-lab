# C-0008 Hybrid Low-Slab Sprint

**Verdict:** `low_slab_beats_terminal_sketch`

**Best (all-IC sketch):** `triad_min_defect_ode` — Ω ≤ **51.0025**
**Best (subclass SOS):** `sos_band123_cubic_ode` — Ω ≤ **17.5363**
(M_sharp=41.2840, M_c0007=41.2840, L-0073 terminal Ω≈54.231)

| Route | Ω(T) worst (low slab) |
|-------|----------------------|
|  **sos_band123_cubic_ode** | 17.5363 |
| sos_greedy_onepol_cubic_ode | 21.7098 |
| triad_min_defect_ode | 51.0025 |
| weighted_min_defect_ode | 57.2967 |
| triad_hybrid_ode | 144.1622 |

Gap vs M_sharp (all-IC): **9.7185**

## SOS cubic (subclass)

- Best band SOS C_ub: 0.8438224893075237
- Best greedy one-pol C_ub: 3.3560651946371634

## Next

1. Prove all-IC C <= C_dagger (fullsym techo blocks naive Shor)
2. Terminal non-partition structure (L-0074 refuted finer clusters)
3. Cluster greedy SOS D=1772 (subclass, ~113h deferred)