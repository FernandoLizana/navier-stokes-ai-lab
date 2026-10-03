# Sprint 02 — Aggressive C-S-0001 attack + L-0006

**Date:** 2026-07-17  
**Route:** B  
**Clay:** none

## Headline

**C-S-0001 numerically refuted (N2, gated).**

| Quantity | Value |
|----------|-------|
| M (C-S-0001) | ≈ 8.895 |
| Best gated ETD | **≈ 32.106** (`abc_n16` shell IC + long polish) |
| Gates | passed |
| Successor | **C-S-0002** with M ≈ **48.159** (= 1.5 × 32.106) |

Shell-supported ABC initial data on N=16, after extended adjoint polish on the **full** grid, exceeds the conjectured enstrophy ceiling.  
This does **not** imply continuum blow-up or refute L-0005 (truncated Galerkin).

## L-0006

Algebraic short-time bound via `‖∇u‖_∞ ≤ √(3M)√(2Ω)`; falls back to L-0001/L-0003.

| Shell | Best closing cap | Winner |
|-------|------------------|--------|
| N12, k≤2 | ≈ **4.43** | (improved vs crude large shells) |
| N12, k≤3 | ≈ 45.1 | still L-0001 (algebraic fails to close at T=0.02) |

## Auditor notes

1. C-S-0001 = shell **IC**, full dealias evolution — high modes allowed after t=0.  
2. L-0005/L-0006 = truncated dynamics — still valid; ABC polish leaves the shell.  
3. Evidence of refutation: N2 (float + integrator/resolution gates), not N9.

## Reproduce

```bash
pytest python/ns_exploration/tests/test_cs0001_attack.py -q
python -m ns_exploration.experiments.sprint02_cs0001_attack
```

## Next

1. Stress **C-S-0002** (M≈48).  
2. Or prove a **truncated-only** corollary: for Galerkin \|k\|≤2, L-0006 already gives Ω≤≈4.43 at these parameters (compare to numerics).  
3. No Clay claims.
