# Sprint 02 Continued Report — NS-MRL

**Date:** 2026-07-16  
**Route:** B only  
**Evidence:** ≤ N2 numerics; conjecture statement N6 (unproved)  
**Clay claim:** none

---

## 1. Implemented

1. **Discrete adjoint** for semi-implicit Euler + analytic VJP of convective nonlinearity (`optimization/adjoint.py`), verified vs FD.
2. **Gate battery** seeds 0–19 (`optimization/candidate_pipeline.py`).
3. **One candidate** after ascent + RK4/ETD + N↔2N cross-check (`datasets/candidates/`).
4. **Conjecture C-0001** — finite-dimensional enstrophy bound only (`conjectures/active/C-0001.json`).

## 2. Tests

```text
pytest python/ns_exploration/tests -q
```

Sprint-2-continued module: 5/5 PASS (VJP↔JVP, instantaneous enstrophy grad, adjoint↔FD, gate smoke, conjecture evaluate).

## 3. Results

| Item | Value |
|------|-------|
| Adjoint vs FD rel error | **0.0018** (agreed) |
| Gate pass rate (20 seeds) | **19/20 = 95%** |
| Candidate seed | 1 (N=12) |
| Integrator RK4↔ETD rel | ~4×10⁻¹⁰ |
| Resolution N↔2N rel | ~6×10⁻⁵ |
| Candidate accepted | **yes** |
| C-0001 status | **exploring** (not refuted) |
| C-0001 probe max enstrophy | ≈ 5.61 ≪ M=20 |

## 4. What this is / is not

- Adjoint is for the **semi-implicit discrete map**, not a continuum NS adjoint theorem.
- Candidate is a **low-resolution exploratory profile**, not a blow-up ansatz.
- C-0001 is **Galerkin-N≤16 only**; Clay implication explicitly **None**.

## 5. Negatives / risks

- 1/20 seeds failed the FD resolution gate — expected; do not optimize those.
- Bound M=20 is loose; easy to leave “exploring” without mathematical progress.
- Full adjoint-based ascent not yet wired into the candidate pipeline (ascent still FD multi-direction).

## 6. Next highest-value step

1. Wire **adjoint gradient ascent** (semi-implicit) with the same N↔2N + integrator gates.
2. Attempt to **refute C-0001** with longer ascent / CMA-ES in coefficient space; if unbroken, **tighten M** to max×1.5 and keep status `exploring`.
3. Only after many failed refutations: attempt a **hand lemma** for the finite Galerkin system (still not Clay).

## Reproduce

```bash
pip install -e "./python[dev]"
pytest python/ns_exploration/tests -q
python -m ns_exploration.experiments.sprint02_continued
```
