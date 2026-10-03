# Sprint — L-0011 Sharpened C-S-0002-SHORT

**Route:** B. **No Clay claim.**

## Rigorous gain (N7)

| | T* for R≈48.16 |
|---|---|
| L-0010 | 0.001011 |
| **L-0011** (div-free + viscous τ) | **0.001240** (~1.23×) |

C-S-0002-SHORT updated to this horizon. N5: `CERT-L0011-…`.

## Conditional L×L majorant (N6 only)

Ignoring H-bilinear remainder: `Ω ≤ B(1+(K U_L τ)²)` gives
`T*≈0.00305` for R≈48 and `Ω(0.02)≈1023` — far below the crude 199650,
but **not** promoted to N7 until cross terms are absorbed.

## Still open

C-S-0002 at `T=0.02`. Bootstrap with high-mode contribution to `||u||_∞`
structurally cannot reach 0.02 for any finite R (`τ W K ≳ 1`).
