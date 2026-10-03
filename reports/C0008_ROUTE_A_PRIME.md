# C-0008 Route A' — Cluster Integral Certificate

Mode: **full**. Clusters: 12. Runtime: 1315.77s.

| Quantity | Value |
|----------|-------|
| I_term^hi (Route A') | 158.0451 |
| Omega(T)^hi | 96.7020 |
| I_*^lo | 1.2993 |
| Per-shell Route A I_hi | 425.3843859528320436290330395846182263747 |
| Ratio A'/A | 0.3715 |
| Closes C-0008 | False |

## Method

Per-cluster min(per-shell sum, B_C max f(r)); combined-output Frobenius for B_C. Still triangle across clusters.

**Conservative per cluster:**
`min(sum_{r in C} B_r f(r), B_C max_{r in C} f(r))`.

Phase E best reference I_hi ~ 78.804 (min(F,1-inf) per shell).