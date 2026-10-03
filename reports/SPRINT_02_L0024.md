# Sprint — L-0024 proves C-0005 (N≤24)

## Result

| claim | status |
|-------|--------|
| Max shell `r=147` has **no** self-triads under dealias | **proved** |
| Spectral defect majorant `(E,Ω)` ODE | **proved** |
| Worst-case `Ω(0.02) ≤ ≈41.74` | **proved** |
| **C-0005** (`M≈61.24`) | **proved** |

Beats Stokes floor only mildly in the majorant; far below envelope 73.5.

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0024.py -q
make sprint2-l0024
```

## Honesty

Finite dealias Galerkin on N≤24 only. Not continuum. Not Clay.
