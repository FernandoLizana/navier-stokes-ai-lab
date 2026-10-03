# L-0076 — N≤24 inclusion sketch (band majorant ladder)

**Status:** `N7_sketch_empirical`  
**C-0008:** `exploring` — does NOT close global claim

## Statement (conditional)

If for each grid `n` the certified terminal majorant `C_term^{hi}(n)` is computed on
the full dealias mode set at that `n` with fixed conventions (2/3 dealias, Hermitian
basis sprintB, terminal weight L-0071), and if mode-set inclusion
`modes(n') ⊆ modes(n)` for `n' < n` implies

\[
\sup_{\|z\|=1} f_n(z) \;\ge\; \sup_{\|z\|\text{ supported on embed}(n')} f_{n'}(z),
\]

then `C_term^{hi}(24)` majorizes any embedded truncation `n' ≤ 24`.

## Verified (empirical, band)

`band_majorant_ladder(radii={1,2,3})` and `{1..6}` check **non-decreasing** sum of per-shell
`C_term_hi` across `n ∈ {12,16,20,24}` under `one_inf_only`.

Implementation: `terminal_weighted/n_coverage.py::band_majorant_ladder`.

**Audit (2026-08-05):** `scripts/run_l0076_n_inclusion_audit.py` →
`reports/L0076_N_INCLUSION_AUDIT.{json,md}`.

| Check | Result |
|-------|--------|
| Wavevector subsets n=12,16,20 ⊆ n=24 | **PASS** |
| Band ladder {1,2,3} monotone | **PASS** (C_term_sum_hi ≈ 16.22, flat) |
| Band ladder {1..6} monotone | **PASS** (C_term_sum_hi ≈ 62.69, flat) |
| Repair full-dealias n=24 Route A | **87/87**, I_hi ≈ **170.59** |

Tests: `ns_exploration/tests/test_l0076_n_inclusion.py` — **3/3 PASS**.

## Still open (N7 → N7)

1. Full-dealias **per-n** ladder at D=6748 for n=12,16,20 (only n=24 complete)
2. Formal embedding of `z_{n'}` into `z_n` preserving terminal functional
3. Proof that mode-set inclusion implies sup monotonicity (ladder is flat empirically, not proved)

## Relation to L-0073R

Cluster residual reweight (L-0073R) applies to temporal integral assembly, independent
of N-domination; both required for rigorous full-dealias certificate.
