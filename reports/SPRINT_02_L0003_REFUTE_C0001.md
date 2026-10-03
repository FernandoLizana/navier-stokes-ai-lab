# Sprint 02 — L-0003 + polish: C-0001 refuted

**Date:** 2026-07-16  
**Route:** B  
**Clay:** none

## Headline

**C-0001 is numerically refuted (N2)** by a CMA candidate polished with adjoint ascent:

| Quantity | Value |
|----------|-------|
| C-0001 M | ≈ 10.994 |
| Gated polished ETD enstrophy | **≈ 13.547** |
| Integrator + N↔2N gates | **passed** |
| Successor | **C-0002** with M ≈ **20.320** (= 1.5 × 13.547) |

This refutes a **finite-dimensional** enstrophy bound only. It does **not** imply continuum blow-up or settle any Clay route.

## L-0003 (proved_restricted, N7)

Uniform-in-time Galerkin bound:

\[
\Omega(t)\le\max\bigl(K^2 E_0,\; 6 M E_0^2/\nu^2\bigr)\qquad\forall t\ge 0.
\]

| N | L-0001 (t=0.02) | L-0002 | L-0003 uniform |
|---|-----------------|--------|----------------|
| 12 | ≈ 1.06×10⁴ | fails to close | ≈ 5.14×10⁴ |
| 16 | ≈ 1.20×10¹¹ | fails | ≫ but finite |

For short T=0.02, L-0001 remains the tighter *closing* ceiling; L-0003 wins for long-time / all-time control at fixed N.

## Files

- `conjectures/proved_restricted/L-0003.json`
- `conjectures/rejected/C-0001.json`
- `conjectures/active/C-0002.json`
- `datasets/candidates/cma_polished_n12.npz`
- `experiments/exploratory/sprint02_l0003_polish/summary.json`

## Evidence ladder

| Claim | Level |
|-------|-------|
| L-0003 | N7 (Galerkin only) |
| C-0001 refutation | N2 (gated float) |
| C-0002 | N6 (proposed) |

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0003_polish.py -q
python -m ns_exploration.experiments.sprint02_l0003_polish
```

## Next

1. Stress-test **C-0002** (M≈20.3) with longer CMA + polish at N=12 and N=16.  
2. Seek embeddings that push L-0001/L-0003 toward O(10)–O(100) for these parameters.  
3. Do not promote N2 refutation to continuum statements.
