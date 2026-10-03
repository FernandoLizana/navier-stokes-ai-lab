# C-0008 Repair Evaluation — Full Dealias one_inf (2026-08-05)

**Conjecture C-0008:** `exploring` — **NOT proved**

## Repair certificate (adversarial PASS)

| Field | Value |
|-------|-------|
| Cert | `certificates/CERT-C0008-repair-full-n24-one-inf.json` |
| Manifest | `python/experiments/terminal_weighted/shell_manifest_repair_full_one_inf_n24.json` |
| Scope | `repair_full` — 87/87 shells, D=6748 |
| Method | `one_inf_only` per-shell MPFR pass |
| Verify | **PASS** |

## Closure vs sharp targets

| Check | Target | Repair bound | Result |
|-------|--------|--------------|--------|
| Integral | I_* ≈ 1.2993 | I_hi ≈ **170.59** | **FAIL** (~131×) |
| Ω(T) | M ≈ 41.284 | ω_T_hi ≈ **101.14** | **FAIL** (~2.4×) |
| `closes_strict` | I_hi < I_* | — | **false** |
| `closes_c0008` | both | — | **false** |

## Historical comparison

| Source | I_hi | Rigorous? |
|--------|------|-----------|
| L-0073 historical (cluster-best) | ~37.92 | **NO** — invalid cluster inequality + non-repair bounds |
| **L-0072 repair (Route A shell sum)** | **~170.59** | **YES** — adversarial verify PASS |

The repair bound is **looser** but **honest**. Historical L-0073 looked better numerically but was not a valid upper bound under repair criteria.

## Lemma supersession

| ID | New status | Repair artifact |
|----|------------|-----------------|
| L-0072 | `superseded_by_repair` | CERT-C0008-repair-full-n24-one-inf.json (Route A) |
| L-0073 | `superseded_historical_invalid` | Route A' cluster **audit PASS** — I_hi ≈ 38.52 (corrected); sketch only |

## Remaining (non-blocking)

- ~~Route A' full cluster manifest + L-0073R cert~~ **DONE** — audit PASS; I_hi ≈ 38.52, does NOT close
- N≤24 inclusion formal proof (L-0076 empirical ladder)
- N7 proof: triangle across clusters + combined B_C pass
- C-0008 remains **exploring** until a bound with I_hi < I_* is proved
