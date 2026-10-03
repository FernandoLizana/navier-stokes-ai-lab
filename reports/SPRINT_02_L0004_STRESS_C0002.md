# Sprint 02 — L-0004 + stress C-0002

**Date:** 2026-07-16  
**Route:** B  
**Clay:** none

## L-0004 (spectral support)

Bounds use `(K_S, M_S)` for fields supported in `|k| ≤ k_shell` instead of the full dealias mask.

| Shell (N=12) | best cap (min L-0001@t, L-0003) |
|--------------|----------------------------------|
| k ≤ 3 | **≈ 45.1** |
| k ≤ 4 | ≈ 425 |

Low-mode families (CMA `k_max≤4`) now sit under **O(10²)** analytic ceilings — closer to numerical targets (~10–20) than full-grid L-0001/L-0003.

Field-conditional `(K_eff, M_eff)` diagnostics tighten further (N2 only).

## C-0002 stress (N=12, N=16; CMA + polish)

| Quantity | Value |
|----------|-------|
| M | ≈ **20.320** |
| Best gated ETD | ≈ **14.267** (`n16_k4`) |
| Refuted? | **No** (margin ≈ **6.05**) |
| Status | `exploring` |

## Evidence

| Item | Level |
|------|-------|
| L-0004 shell bounds | N7 (hypothesis H(S)) |
| C-0002 survival | N2 |
| C-0002 statement | N6 unproved |

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0004_stress_c0002.py -q
python -m ns_exploration.experiments.sprint02_stress_c0002
```

## Next

1. Formulate **C-0002-shell**: enstrophy bound for `|k|≤4` IC families only — may be provable with L-0004 + sharper embedding.
2. Push polish/CMA at N=16 with more generations before tightening M below ~15.
3. Still no continuum / Clay claims.
