# Sprint — L-0021 mono-radial quartic → C-R-0001

## Result

| claim | status |
|-------|--------|
| Quartic Gram `α_★ = max_r α_r = 14` (r=74) | **proved** (L-0021) |
| Mono-radial Duhamel `Ω(0.02) ≤ ≈71.73` | **proved** (L-0020+21) |
| **C-R-0001** (single radial shell IC) | **proved** |
| All-IC **C-0005** (M≈61.24) | still **open** (`N_crit≈9.67`; `α_★=14` too big) |

Beats envelope 73.5 on the mono-radial subclass.

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0021.py -q
make sprint2-l0021
```

## Honesty

Finite Galerkin, mono-radial ICs only. Quartic Gram is a Frobenius relaxation
(upper-bounds the true shell max). Not all-IC. Not Clay.
