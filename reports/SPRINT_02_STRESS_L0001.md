# Sprint 02 Stress C-0001 + L-0001 Report

**Date:** 2026-07-16  
**Route:** B  
**Clay:** none

## What was done

1. Longer adjoint stress of C-0001 at **N=12** and **N=16**.
2. Hand lemma **L-0001**: explicit exponential enstrophy ceiling for Fourier–Galerkin (N7, restricted).

## Stress results

| Quantity | Value |
|----------|-------|
| M tested | ≈ 10.994 |
| Best ETD enstrophy found | ≈ **8.491** |
| Refuted with gates? | **No** |
| C-0001 status | `exploring` (M unchanged) |

## L-0001 (proved_restricted)

\[
\Omega(t)\le K^2 E_0\exp\bigl(2K\sqrt{3M}\sqrt{2E_0}\,t\bigr)
\]

| N | K | M | Ω(0.02) cap |
|---|-----|------|-------------|
| 12 | (from run) | | ≈ 1.06×10⁴ |
| 16 | ≈8.66 | 1331 | ≈ 1.20×10¹¹ |

Far above C-0001’s numerical target (~11). Useful as a **constructive finite-N ceiling**, useless for Clay (blow-up of constants in N).

## Evidence

- Stress campaign: **N2**
- L-0001: **N7** (elementary estimates on the Galerkin ODE)
- C-0001 explicit M: still **N6** / unproved

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0001_and_stress.py -q
python -m ns_exploration.experiments.sprint02_stress_c0001
```

## Next highest-value step

1. Improve L-0001 using viscosity / Gevrey / better embeddings to shrink the gap toward C-0001.
2. Or attack C-0001 with coefficient-space CMA-ES under energy constraint.
3. Do **not** claim continuum consequences.
