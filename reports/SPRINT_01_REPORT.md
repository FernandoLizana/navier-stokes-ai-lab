# Sprint 01 Report — NS-MRL

**Date:** 2026-07-16  
**Route:** B (primary). Secondary D not executed.  
**Claim level:** No Clay resolution. No N9+ proofs.

---

## 1. What was implemented

- Full repository skeleton (docs, literature, python, julia, lean stubs, experiments, datasets, …).
- Official Clay target documentation (`docs/OFFICIAL_CLAY_TARGET.md`).
- Mathematical foundations, function spaces, regularity criteria docs (classical citations only).
- Literature matrix with **64** primary/classified rows (`literature/literature_matrix.csv`).
- Python exploratory pseudospectral solver on \(\mathbb{T}^3=[0,2\pi)^3\):
  - unified Fourier conventions;
  - Leray projector;
  - 2/3 dealiasing;
  - integrators RK4, ETD-RK2, semi-implicit Euler;
  - diagnostics (energy, enstrophy, helicity, spectral tail, analyticity proxy);
  - checkpoints (npz).
- Initial conditions: Taylor–Green, ABC, random div-free, periodic vortex tubes.
- Finite-mode interval arithmetic module (Python); Julia `ValidatedNS` skeleton (Julia not installed on this machine).
- Lean 4 file stubs (no theorems).
- Pilot campaign (infrastructure only).

## 2. Tests that pass

All **15** critical tests in `python/ns_exploration/tests/test_sprint01.py`:

| Test | Result |
|------|--------|
| FFT roundtrip | PASS |
| Leray kills divergence | PASS |
| TG / ABC / tubes div-free | PASS |
| Linear diffusion exact mode | PASS |
| Linear diffusion via integrators | PASS |
| Energy decays (f=0) | PASS |
| Divergence preserved in evolution | PASS |
| Nonlinear energy channel (ABC) | PASS |
| Mesh convergence (linear) | PASS |
| Integrator agreement short time | PASS |
| Interval Leray / modes / convolution / heat | PASS |

```text
pytest python/ns_exploration/tests/test_sprint01.py -v
# 15 passed
```

## 3. Errors found and fixed

1. **Vortex-tube IC:** `enforce_reality` reintroduced \(O(10^{-10})\) divergence. Fixed by re-applying Leray after the reality round-trip.
2. **Empty repo:** no legacy code conflicts; greenfield build.
3. **Julia:** not available in PATH; Julia interval tests not executed here. Python interval tests cover the Sprint 1 finite-mode probes.

## 4. Results that are floating-point (exploratory)

- All pseudospectral evolutions (pilot + solver tests beyond exact linear mode algebra).
- Diagnostics (enstrophy, \(\|\omega\|_\infty\), analyticity proxy).
- Pilot campaign (72 runs): evidence **N2**.

Pilot summary:

```json
{
  "n_runs": 72,
  "max_div": 9.71e-17,
  "n_stopped": 0,
  "purpose": "infrastructure_validation",
  "not_a_singularity_search": true
}
```

## 5. Results that are validated (narrow sense)

- **Finite-mode** interval checks (Python): Leray divergence enclosure, short convolution enclosure, linear heat factor enclosure.
- Evidence target **N5 for those finite identities only**.
- **Not** a validated continuum NS solution; **not** Fourier-tail control; **not** N9.

## 6. What has **not** been proved

- Global regularity on \(\mathbb{T}^3\) (Route B).
- Finite-time singularity (Route D/C).
- Any new analytical inequality.
- Continuum a posteriori existence certificates.
- Lean theorems.
- That numerical spectral decay implies smoothness for all time.

## 7. Priority mathematical risk

**Confusing successful float evolution / small divergence with continuum regularity**, and/or later confusing Galerkin blow-up with PDE blow-up. Mitigation: evidence ladder, stop conditions, explicit Route tags, auditor checklist in `docs/SCIENTIFIC_INTEGRITY.md`.

Secondary risk: **un-audited interval FFT** if introduced prematurely — keep exploratory FFT separate from validated path.

## 8. Highest-value next experiment

**Do not start full adversarial enstrophy maximization yet** beyond a controlled dry-run until:

1. Add manufactured / Beltrami quantitative residual report at two resolutions with documented dealiasing error;
2. Implement an independent **2D embedded** or low-mode Galerkin cross-check of energy balance;
3. Then Sprint 2 item 1–2: adjoint/finite-difference gradient of enstrophy at fixed energy on **low N**, with mandatory resolution doubling before any conjecture.

---

## Answers to mandatory Sprint questions (compressed)

1. **Implemented:** docs + literature + exploratory spectral NS + ICs + interval probes + pilot.  
2. **Tests:** 15/15 PASS.  
3. **Errors:** reality→div leak on tubes (fixed); Julia missing.  
4. **Floating:** solver + pilot N2.  
5. **Validated:** finite-mode interval identities only (narrow N5).  
6. **Not demonstrated:** Clay A–D.  
7. **Risk:** numerics ≠ proof; truncation ≠ continuum.  
8. **Next:** residual/cross-check, then cautious adjoint enstrophy search with resolution gates.

## Reproduce

```bash
pip install -e "./python[dev]"
pytest python/ns_exploration/tests/test_sprint01.py -v
python -m ns_exploration.experiments.pilot_campaign
# optional when Julia is installed:
# make interval-test
```

## Evidence levels used this sprint

| Artifact | Level |
|----------|-------|
| Documentation of classical theorems | N7 (citations) |
| Solver / pilot | N2 |
| Interval finite-mode probes | N5 (finite only) |
| Clay solution claim | **none / forbidden** |
