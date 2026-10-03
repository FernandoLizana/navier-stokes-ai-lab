# Sprint 02 Partial Report — NS-MRL

**Date:** 2026-07-16  
**Route:** B only  
**Scope:** Post–Sprint-1 gates before adversarial search  
**Claim:** No Clay resolution. Evidence ≤ N2 (N3 pending independent reproduction).

---

## 1. What was implemented

1. **Residual suite** (`diagnostics/residuals.py`)
   - ABC Beltrami nonlinear residual at N=16,32 with/without dealiasing
   - Manufactured single-mode viscous decay residual at N=16,32
2. **Energy cross-checks** (`diagnostics/energy_crosscheck.py`)
   - Parseval spectral vs physical
   - Independent 2D-reduced nonlinear operator vs full 3D on z-invariant TG
   - Stepwise energy vs dissipation proxy
3. **Enstrophy FD gradient + resolution gate** (`optimization/enstrophy_gradient.py`)
   - Fixed-energy sphere projection
   - Central FD directional derivative
   - Mandatory N vs 2N agreement gate
   - Single ascent step only if gate passes
4. Experiment runner + tests + this report

## 2. Tests

```text
pytest python/ns_exploration/tests -q
# 23 passed
```

Sprint-2-specific: 8/8 PASS in `test_sprint02_partial.py`.

## 3. Quantitative results (from `experiments/exploratory/sprint02_partial/summary.json`)

| Check | Result |
|-------|--------|
| Beltrami residual rel (N=16, dealias) | ~7×10⁻¹⁶ |
| Beltrami residual rel (N=32, dealias) | ~1×10⁻¹⁵ |
| Manufactured mode rel error | ~7×10⁻¹⁶ |
| Parseval rel diff | 0 |
| 2D vs 3D operator rel diff | 0 |
| Energy step residual (max abs) | ~5×10⁻⁹ |
| FD gate relative agreement (N=12 vs 24) | **0.015** (tol 0.35) → **PASSED** |
| One ascent step ΔJ | +0.029 (N=12 only) |

## 4. Floating vs validated

- All of the above are **floating-point exploratory (N2)**.
- No new interval certificates in this partial sprint.
- Gate pass does **not** authorize conjectures about maximal enstrophy or blow-up.

## 5. Errors / negatives

- None blocking. Gate passed on the chosen seed/time window.
- Adjoint (exact reverse-mode) not implemented yet — FD only.
- Ascent at N=12 is **not** resolution-validated as an optimizer trajectory (only the directional derivative gate was).

## 6. Mathematical risks (updated)

1. Optimizing enstrophy on under-resolved grids (mitigated by gate; still weak).
2. Interpreting ΔJ>0 as approach to singularity — **forbidden**.
3. 2D embedded agreement does not lift to 3D vortex stretching control.

## 7. What was demonstrated / not

**Demonstrated (numéricamente):** residual infrastructure healthy; Parseval OK; 2D cross-operator OK; short-time FD enstrophy derivative agrees at N and 2N within 2% for one seed.

**Not demonstrated:** maximizers, blow-up, regularity, adjoint correctness, continuum certificates.

## 8. Next highest-value step

1. Implement **adjoint** enstrophy gradient (or verify FD vs adjoint on low modes).
2. Multi-seed gate battery (seeds 0–19) — proceed to ascent campaigns only on seeds that pass N↔2N.
3. Extract **one** candidate profile only after gate + integrator cross-check (RK4 vs ETD).
4. Formulate a **refutable** finite-dimensional conjecture (e.g. bound on enstrophy growth for Galerkin-N≤N0) — not Clay.

## Reproduce

```bash
pip install -e "./python[dev]"
pytest python/ns_exploration/tests -q
python -m ns_exploration.experiments.sprint02_partial
```
