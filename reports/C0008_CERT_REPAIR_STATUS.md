# C-0008 Certification Repair — Frozen Baseline

**Date:** 2026-08-05  
**Conjecture C-0008:** `exploring` — **NOT proved**  
**Repair phase:** 12 — **full-dealias COMPLETE + cert PASS**

## Phase 12 — full-dealias closure (2026-08-05)

| Item | Status |
|------|--------|
| 87/87 shells one_inf_only | **DONE** |
| `shell_manifest_repair_full_one_inf_n24.json` | **DONE** |
| `CERT-C0008-repair-full-n24-one-inf.json` | **DONE** — verify **PASS** |
| Evaluation `reports/C0008_REPAIR_EVALUATION.md` | **DONE** |
| L-0072 / L-0073 registry supersession | **DONE** |
| Route A' cluster manifest L-0073R | **DONE** — `cluster_manifest_repair_full_n24_one_inf.json` (11 clusters) |
| Cluster cert L-0073R | **DONE** — verify PASS; audit PASS (`L0073R_CLUSTER_AUDIT.md`); I_hi ≈ **38.52** |
| L-0076 N≤24 inclusion audit | **DONE** — `reports/L0076_N_INCLUSION_AUDIT.md` |
| Band {1..6} rational regen | **DONE** — `shell_manifest_repair_band6_one_inf_rational.json` (~65s); shell blocks match legacy within MPFR noise |
| Full n=24 rational regen | **RUNNING** — `regenerate_c0008_full_dealias_one_inf.py --rational --watchdog` (87 shells, days) |

### Repair bounds (rigorous, adversarial PASS)

| | Repair | Target | Closes? |
|---|--------|--------|---------|
| I_hi | **170.59** | I_* ≈ 1.2993 | **NO** (~131×) |
| ω_T_hi | **101.14** | M ≈ 41.284 | **NO** (~2.4×) |

### Cluster Route A' (post-audit 2026-08-06)

| | Cluster L-0073R (corrected) | Route A | Target | Closes? |
|---|---------------------------|---------|--------|---------|
| I_hi | **38.52** | 170.59 | I_* ≈ 1.299 | **NO** (~30×) |
| ω_T_hi | ~54.2 | 101.14 | M ≈ 41.28 | **NO** (~1.3×) |

Bug fixed: `cluster_residual_contrib_hi` erroneously multiplied residual by `max_phi` (~50× deflation). Corrected bound ≈ historical ~37.92 but still **not** rigorous N7 closure.

## Historical bounds (superseded / invalid)

| ID | Previous label | Repair status |
|----|----------------|---------------|
| L-0072 | route_refuted | **`superseded_by_repair`** — CERT repair full one_inf PASS |
| L-0073 | best_known_bound | **`superseded_historical_invalid`** — repair Route A PASS; historical ~37.92 NOT rigorous |
| L-0074 | partition_ladder_refuted | **`evaluated_configurations_did_not_close`** |
| L-0075 | best_low_slab_hybrid | **`N2/N6_sketch`** — grid ODE, not N7 |

## P0 defects identified

1. **Frobenius streaming** sums `|delta|^2` per update, not `|M_ij|^2` per entry — **FIXED** (entry aggregation then square; re-enabled for band)
2. **Cluster max(φ)·B_C** inequality invalid for matrix sums — **REPLACED** by residual reweight (L-0073R)
3. **Verifiers** trusted declared totals; missing manifest fail-soft — **FIXED**
4. **Intervals** used float constants and wrong sqrt orientation for Omega — **FIXED**
5. **direct-vs-tensor test** compared `f_from_G(G,z)` vs itself — **FIXED** (4 paths)
6. **Permutation set** in streaming M collapses multiplicities — **FIXED** (6-tuple list in tensor + l0045)
7. **Coefficient construction** float-first — **MPFR-directed G_lo/G_hi** encloses float G; **`rational_basis.py`** verifies pol/Leray vs float reference (worst rel ≲ 10⁻¹⁶)
8. **Adaptive temporal** early-stop before queue empty — **FIXED**
9. **L-0075** mislabeled N7; M_c0007 conflated with M_target — **FIXED**
10. **N≤24 coverage** remains `N7_sketch` — formal proof deferred

## Active bound method (repair baseline)

```
bound_method = one_inf_only
evidence_level = N4 (pending exact coefficient enclosure)
rational_coefficients = **production on band {1..6}** (`shell_manifest_repair_band6_one_inf_rational.json`); full n=24 rational regen **pending** (days, external shell)
```

Do **not** run full-dealias N=24 until GO checklist in `C0008_CLUSTER_REWEIGHT_REPAIR.md` passes.

## Phase 5–6 progress (2026-08-04)

| Item | Status |
|------|--------|
| `rational_basis.py` (MPFR pol/Leray from integer k) | **DONE** — 6 modes, worst rel ≲ 10⁻¹⁶ |
| `certified_cluster_one_inf_hi` (combined `‖Σ M_r‖` pass) | **DONE** |
| `cluster_residual_prototype` uses combined B_C | **DONE** |
| Band `{1,2,3,4,5,6}` four-path + interval + cluster | **DONE** (smoke) |
| Tests | **18/18 PASS** |

### Cluster residual (one_inf, n=24)

| Band | B_C combined | B_C Σ B_r (triangle) | best residual α=min |
|------|-------------|----------------------|---------------------|
| `{1,2,3}` | 13.21 | 16.22 | 0.00525 |
| `{1..6}` | 32.59 | 62.69 | 0.01295 |

Combined pass is strictly tighter than triangle sum (as expected). Legacy `B_C×max(φ)` remains invalid (~0.26 / ~0.65).

### Remaining NO-GO blockers

- ~~Full dealias N=24 manifest/certificate at D=6748~~ **DONE**
- N≤24 inclusion formal proof (L-0076 band ladder **empirical PASS**)
- N7: triangle across clusters + exact coefficient enclosure
- C-0008 sharp closure (I_hi < I_*): **OPEN** — best repair cluster ~38.5, Route A ~170.6

## Phase 7 progress (2026-08-04)

| Item | Status |
|------|--------|
| Rational interval tensor (`mpfr_leray` + `mpfr_vdot_real`) | **DONE** — midpoint gap ≲ 10⁻¹⁶ vs float interval |
| Frobenius entry aggregation repair | **DONE** — re-enabled on band |
| `shell_manifest_repair_band6_one_inf.json` | **DONE** |
| Tests repair+band | **20/20 PASS** |

## Phase 9–10 progress (2026-08-04)

| Item | Status |
|------|--------|
| `repair_full_manifest.py` — per-n checkpoints + watchdog | **DONE** |
| Legacy `upgrade_shell_ck` migration for resume | **DONE** (shell 1 @ 34M/45.5M pairs) |
| `build_cluster_manifest` L-0073R for `one_inf_only` | **DONE** |
| `build_c0008_repair_full_certificate.py` | **DONE** |
| n=24 watchdog | **DONE** (87/87 shells) |
| Tests | **30+ PASS** |

## Phase 11 — superseded by Phase 12

Historical notes below retained for audit trail only.

## Phase 11 progress (2026-08-04)

| Item | Status |
|------|--------|
| `c0008_repair_pipeline.py` — enrich → cluster → cert → verify | **DONE** |
| `scripts/finalize_c0008_repair.py` | **DONE** |
| `scripts/bootstrap_c0008_repair.ps1` + `run_c0008_repair_all.ps1` | **DONE** |
| `requirements-c0008-repair.txt` | **DONE** |
| CI `.github/workflows/c0008-repair-band.yml` | **DONE** |
| Band cert adversarial verify | **PASS** (`verify_ok=True`) |
| n=24 shell 1 | **DONE** (C_term_hi≈9.02); shell 2 **IN PROGRESS** |
| Full-dealias GO | **NO-GO** until 87/87 shells + `repair_full` cert PASS |

### Pipeline note

Declared `integral_interval` stays **Route A** (shell-sum, verifier-recomputable). Route A' cluster sketch stored in `best_route_sketch` only — does **not** supersede declared bound until cluster manifest verification is wired.

### Monitor full-dealias regeneration (adaptive max-CPU mode)

```powershell
Set-Location python
python scripts/regenerate_c0008_full_dealias_one_inf.py --n 24 --status
Get-Content experiments\terminal_weighted\repair_one_inf_heartbeat_n24.txt
Get-Content experiments\terminal_weighted\repair_watchdog_state_n24.json
```

Resilient launcher (8 workers adaptive, auto-finalize):

```powershell
powershell -File scripts\run_c0008_repair_full_resilient.ps1 -Workers 8
```

Checkpoints: per-shell `shell_one_inf_ck/n24/shell_NNN.json` + manifest `repair_one_inf_checkpoint_n24.json`. Watchdog backs off workers on OOM/crash.
