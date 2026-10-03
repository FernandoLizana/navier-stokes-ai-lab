# Sprint — L-0014 Triad multiplicity → longer C-S-0002-SHORT

## Result

| bound | T* for R≈48.16 (N=16, ℓ∞≤4) |
|---|---|
| L-0013 (vector CS) | ≈0.003235 |
| **L-0014** (shell→high `R_H=324`) | **≈0.00437** (~1.35×; ~22% of T=0.02) |

C-S-0002-SHORT updated. N5: `CERT-L0014-triad-CS0002-short-N16`.

## Idea

`‖P_H N‖_ℓ₂ ≤ √R_H ‖û‖₂ ‖(|q|û)‖₂ = 2√(R_H E Ω)` replaces
`‖N‖ ≤ ‖u‖_∞‖∇u‖₂` (which costs `√M_L`). Effective `U_eff=√R_H√(2E0)=18`.

## Commands

```
pytest python/ns_exploration/tests/test_l0014.py -q
make sprint2-l0014
```

## Honesty

Finite Galerkin, short horizon only. Not continuum. Not a Clay claim.
Parent C-S-0002 @ T=0.02 still open.
