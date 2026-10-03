# L-0006 — Improved short-time shell bound

**Evidence:** N7 (restricted)  
**Clay:** none

## Idea

Replace `||∇u||_∞ ≤ K√(3M)√(2E)` by `||∇u||_∞ ≤ √(3M)√(2Ω)` and integrate
`dΩ/dt ≤ 2√2√(3M) Ω^{3/2}` to get an algebraic short-time ceiling when it closes;
otherwise fall back to L-0001 / L-0003 and take the best finite bound.

Still shell Galerkin only — not C-S-0001 (full-grid evolution).
