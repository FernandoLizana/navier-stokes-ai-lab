# Sprint — L-0016 Radial Young → longer C-S-0002-SHORT

## Result

| bound | T* for R≈48.16 (N=16, ℓ∞≤4) |
|---|---|
| L-0015 (`S_max`) | ≈0.005967 |
| **L-0016** (`ρ_★=1392`) | **≈0.00912** (~1.53×; ~46% of T=0.02) |

C-S-0002-SHORT updated. N5: `CERT-L0016-radial-CS0002-short-N16`.

## Idea

Young `‖A∗B‖₂≤‖A‖₂‖B‖₁` + radial shells: `N_max=2 E0 √(max_r r m_r)`.
Worst shell `r=29` with `m_r=48` ⇒ `U_eff≈5.39`.

## Commands

```
pytest python/ns_exploration/tests/test_l0016.py -q
make sprint2-l0016
```

## Honesty

Finite Galerkin, short horizon only. Not continuum. Not a Clay claim.
Parent C-S-0002 @ T=0.02 still open (H×L `W` term now dominates).
