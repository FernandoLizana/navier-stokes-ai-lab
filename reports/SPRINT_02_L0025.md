# Sprint — L-0025 proves C-0006 (N≤32); opens C-0007

## Result

| claim | status |
|-------|--------|
| L-0024 method at N=32 (no max-shell self-triads) | **works** |
| **C-0006** (N≤32, M=1.5×Stokes_floor≈67.77) | **proved** (majorant ≈46.15) |
| **C-0007** (N≤24, M≈41.28 in floor–majorant gap) | **opened** (active) |

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0025.py -q
make sprint2-l0025
```

## Honesty

Finite dealias Galerkin only. Not continuum. Not Clay.
