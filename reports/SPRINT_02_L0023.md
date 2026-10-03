# Sprint — L-0023 cubic dissipation (new path to C-0005)

## Result

| claim | status |
|-------|--------|
| ODE: `z' ≤ -(ν/E0) z³ + (C/2) z²` from `Σ|k|⁴ ≥ 2Ω²/E` | **proved** (framework) |
| `C ≤ C_★ ≈ 2.154` ⇒ `Ω(0.02) ≤ M_{C-0005}` from `Ω0 ≤ 73.5` | **proved** (arithmetic) |
| Empirics `C_emp ≪ C_★` | **N2 support** |
| Proved `C ≤ C_★` | **open** |
| **C-0005** | still **open** |

Uniform Duhamel (L-0020/22) is a dead end for all-IC. This ODE path has
~50× slack vs observed stretch constants.

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0023.py -q
make sprint2-l0023
```

## Honesty

Finite Galerkin only. Empirics ≤ N2. Not Clay.
