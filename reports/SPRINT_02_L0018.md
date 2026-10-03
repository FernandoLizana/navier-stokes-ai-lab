# Sprint — L-0018 Spectral envelope closes C-S-0002

## Result

| bound | claim |
|---|---|
| L-0017 (cancellation ODE, loose Ω bookkeeping) | T*≈0.01606; parent @0.02 open under that ODE |
| **L-0018** (Parseval envelope `Ω ≤ K² E`) | **`Ω(t) ≤ 37.5 < M≈48.16` for all `t≥0` on N≤16** |

**C-S-0002 proved** (restricted: finite dealiased Galerkin, N≤16). N5: `CERT-L0018-envelope-CS0002-N16`.

## Idea

On the dealias mask, `Ω = (1/2) Σ |k|²‖û‖² ≤ K²(N) E` with `K²(16)=75`. Energy decay `E(t)≤E0` ⇒ `Ω≤37.5`. Shell IC is a subclass.

The L-0012–17 comparison used `Ω ≤ K_IC² E0 + K_full² E_H` (double energy budget); that artifact created a false gap to T=0.02.

## Does not close

- **C-0002** (`M≈20.32`): envelope 37.5 is too weak.
- Continuum / Clay: no.

## Commands

```
pytest python/ns_exploration/tests/test_l0018.py -q
make sprint2-l0018
```
