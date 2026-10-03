# Sprint 02 Adjoint Ascent Report — NS-MRL

**Date:** 2026-07-16  
**Route:** B  
**Evidence:** N2 numerics; C-0001 statement N6 (explicit M still unproved)  
**Clay:** none

## Implemented

1. `optimization/adjoint_ascent.py` — energy-sphere ascent via discrete adjoint + backtracking; RK4/ETD and N↔2N gates.
2. `conjectures/c0001_adversarial.py` — adversarial refutation campaign; M tightening; compactness lemma note.
3. Experiment `sprint02_adjoint_ascent` + tests.

## Results

| Item | Value |
|------|-------|
| Demo ascent (seed 1) | J: 5.14 → **7.28** |
| Demo gates | **passed** |
| Campaign best ETD enstrophy | **≈ 7.33** |
| Refuted M=20 with gates? | **No** |
| New explicit M | **≈ 10.994** (= 1.5 × best) |
| C-0001 status | `exploring` |

## Mathematical note (finite-N)

For each fixed Galerkin resolution N, the energy sphere of divergence-free fields is compact in finite dimensions and the discrete flow map is continuous, so \(\max \mathcal{E}(T)\) exists. This **trivial lemma** does not give an explicit M, uniformity in N, or any Clay statement. C-0001 remains an **explicit-constant** finite-dimensional conjecture.

## Tests

```bash
pytest python/ns_exploration/tests/test_sprint02_adjoint_ascent.py -q
python -m ns_exploration.experiments.sprint02_adjoint_ascent
```

## Next

1. Longer adjoint ascent / multi-start to stress-test new M≈11.
2. Repeat campaign at N=16 (still finite).
3. If M survives: attempt a **hand upper bound** for the Galerkin system using energy + enstrophy ODEs with explicit constants (still not Clay).
