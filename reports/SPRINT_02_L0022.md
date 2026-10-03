# Sprint — L-0022 gamma bootstrap (C-0005 path closed as uniform)

## Result

| claim | status |
|-------|--------|
| Arithmetic: C-0005 closes if `γ ≤ γ_crit≈1.7469` | **proved** (conditional) |
| Uniform `‖N‖≤N_crit` / `γ≤γ_crit` for all-IC | **false** (adversarial N2) |
| **C-0005** refuted? | **no** (`Ω(T)_emp ≲ 41 ≪ 61`) |
| **C-0005** proved? | still **open** |

Adversarial VJP ascent reaches `‖N‖ ≳ 35–40 ≫ N_crit≈9.67` and
`γ ≳ 10 ≫ γ_crit`, so L-0020 with a **uniform** `N_*` cannot close all-IC
C-0005. Those same fields keep `Ω(0.02) ≪ M` — Duhamel is pessimistic.

Prior “`‖N‖≲2`” applied only to unstructured random ICs, not adversaries.

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0022.py -q
make sprint2-l0022
```

## Next

Need a sharper-than-uniform majorant (structure/cancellation in `N`, or a
true enstrophy ODE with better production control). Mono-radial C-R-0001
remains the best proved gap-below-envelope result.

## Honesty

Finite Galerkin only. Empirics ≤ N2. Not Clay.
