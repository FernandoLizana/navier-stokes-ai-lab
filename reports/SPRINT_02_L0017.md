# Sprint — L-0017 Radial H×L cross → longer C-S-0002-SHORT

## Result

| bound | T* for R≈48.16 (N=16, ℓ∞≤4) |
|---|---|
| L-0016 (radial L×L only) | ≈0.009116 |
| **L-0017** (radial L×L + H×L) | **≈0.01606** (~1.76×; **~80%** of T=0.02) |

C-S-0002-SHORT updated. N5: `CERT-L0017-cross-CS0002-short-N16`.

## Idea

Young on `(u_H·∇)u_L`: cross `z`-coeff `γ_cross=√(2ρ_★ E0)≈37.3` replaces `W√B≈170`.

## Commands

```
pytest python/ns_exploration/tests/test_l0017.py -q
make sprint2-l0017
```

## Honesty

Finite Galerkin, short horizon only. Not continuum. Not a Clay claim.
Parent C-S-0002 @ T=0.02 still open (`Ω(0.02)≈67.6 > R≈48.2`).
