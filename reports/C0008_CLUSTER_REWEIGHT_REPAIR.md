# C-0008 Cluster Reweight Repair

**Date:** 2026-08-04  
**C-0008:** `exploring`

## Retired cluster inequality

`cluster_integral_hi = min(per_shell_sum, B_C * max_factor)` used invalid
`max(phi)*B_C` for combined weighted operator.

## Replacement (L-0073R)

See `lemmas/L-0073R-cluster-residual-reweight.md`.

Implementation: `cluster_bounds.cluster_residual_upper()`.

## Status changes

| Lemma | Old | New |
|-------|-----|-----|
| L-0073 | best_known_bound | **invalid_pending_cluster_reweight_proof** |
| L-0074 | partition_ladder_refuted | **evaluated_configurations_did_not_close** |

Historical manifests preserved; not deleted.

## GO checklist before N=24 full-dealias

- [x] Verifiers adversarial (declared < recomputed fails)
- [x] Frobenius disabled; one_inf_only baseline
- [x] Interval sqrt orientation fixed
- [x] Four-path direct-vs-tensor test
- [x] Permutation multiplicity fix in streaming M
- [x] Adaptive temporal no early-stop
- [x] Exact {1..3} coefficient enclosure prototype (rational Leray + MPFR vdot)
- [x] Cluster residual B_C combined pass (one_inf)
- [ ] N≤24 inclusion proof (N7_sketch → N7) — **L-0076** band ladder empirical OK

**Decision:** **DONE** for n=24 full-dealias (Phase 12). Remaining: N7 coefficient enclosure + cross-cluster proof.
