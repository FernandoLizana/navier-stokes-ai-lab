# Sprint 02 — L-0002 + CMA-lite vs C-0001

**Date:** 2026-07-16  
**Route:** B  
**Clay:** none

## Implemented

1. **L-0002** — viscous Galerkin comparison ODE with \(\|\nabla u\|_\infty\le\sqrt{3M}\sqrt{2\Omega}\) and \(\|\nabla\omega\|^2\ge 2\Omega\).
2. **CMA-lite** coefficient-space \((\mu/\lambda)\)-ES attacking C-0001 on low modes.

## Results

| Item | Value |
|------|-------|
| CMA best ETD enstrophy | **≈ 10.236** |
| C-0001 M | ≈ 10.994 |
| Refuted? | **No** (margin ≈ 0.76) |
| L-0002 at N=8,12,16 | **estimate fails to close** (\(t_* < 0.02\)) |
| L-0001 N=12 ceiling | ≈ 1.06×10⁴ (still valid, loose) |

### Auditor note (critical)

L-0002’s comparison ODE blowing up is **failure of the estimate**, not evidence of Galerkin or PDE singularity. Documented explicitly in lemma notes.

CMA closed much of the gap to C-0001 (8.49 → 10.24) without crossing M.

## Evidence

| Artifact | Level |
|----------|-------|
| L-0002 | N7 (restricted; often non-closing) |
| CMA attack | N2 |
| C-0001 | N6, still unproved |

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0002_cmaes.py -q
python -m ns_exploration.experiments.sprint02_l0002_cmaes
```

## Next highest-value step

1. Tighten embeddings for L-0002/L-0003 (e.g. restrict \(M\) to active energy-containing shells, or Gevrey weights) so a **closing** bound sits nearer C-0001.
2. Or push CMA with higher \(k_{\max}\) + resolution gate before claiming near-maximizer.
3. If numerical max stabilizes near ~10.3, reformulate C-0001 with \(M=1.1\times\) that value as a sharper conjecture still finite-dimensional.
