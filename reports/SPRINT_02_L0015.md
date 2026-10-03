# Sprint — L-0015 S_max → longer C-S-0002-SHORT

## Result

| bound | T* for R≈48.16 (N=16, ℓ∞≤4) |
|---|---|
| L-0014 (`R_H=324`) | ≈0.004367 |
| **L-0015** (`S_max=6132`) | **≈0.00597** (~1.37×; ~30% of T=0.02) |

C-S-0002-SHORT updated. N5: `CERT-L0015-smax-CS0002-short-N16`.

## Idea

`‖P_H N‖² ≤ 4 E² S_max` with `S_max = max_p Σ_{q:p+q∈H} |q|²`,
so `‖P_H N‖ ≤ 2 E0 √S_max` (energy-only; no Ω). Effective `U_eff≈11.30`.

## Commands

```
pytest python/ns_exploration/tests/test_l0015.py -q
make sprint2-l0015
```

## Honesty

Finite Galerkin, short horizon only. Not continuum. Not a Clay claim.
Parent C-S-0002 @ T=0.02 still open (need ~2× sharper `U_eff`).
