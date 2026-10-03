# L-0004 — Spectral support restricted Galerkin bounds

**Evidence:** N7 (shell hypothesis); field-conditional diagnostics N2  
**Clay:** none

## Idea

Replace global `(K,M)` from full dealias mask with `(K_S, M_S)` for fields supported in `|k|≤k_shell`.
CMA/ES IC families with `k_max≤4` fall under these hypotheses.

## Example caps (E0=0.5, t=0.02, ν=0.1) — see experiment summary for full table

Low shells can yield **O(10)–O(10³)** ceilings vs **O(10⁴)+** for full grid.

Field-conditional `(K_eff, M_eff)` can be tighter still (diagnostic only).
