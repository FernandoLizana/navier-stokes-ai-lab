# Sprint — L-0013 Vector CS embedding → longer C-S-0002-SHORT

## Result

| bound | T* for R≈48.16 (N=16, ℓ∞≤4) |
|---|---|
| L-0012 (NS cancellations) | ≈0.002284 |
| **L-0013** (vector CS `‖u‖_∞≤√M√(2E)` + L-0012 ODE) | **≈0.00324** (~1.42×) |

C-S-0002-SHORT updated. N5: `CERT-L0013-embed-CS0002-short-N16`.

Parent **C-S-0002 @ T=0.02** still open. Even with W=0, need `U_L≲5.94` vs current ≈26.98.

## Fix

L-0011 used `√(dof·M)√(2E)` with dof=2, i.e. an extra √2 over the vector mode CS
`Σ‖û_k‖ ≤ √M √(Σ‖û_k‖²) = √M √(2E)`. Also exclude `k=0` from `M_L` (728 vs 729).

## Commands

```
pytest python/ns_exploration/tests/test_l0013.py -q
make sprint2-l0013
```

## Honesty

Finite Galerkin, short horizon only. Not continuum. Not a Clay claim.
