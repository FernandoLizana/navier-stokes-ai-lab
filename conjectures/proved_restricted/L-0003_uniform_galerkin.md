# L-0003 — Uniform-in-time Galerkin enstrophy bound

**Status:** proved_restricted  
**Evidence:** N7  
**Clay:** none

## Statement

For mean-zero Fourier–Galerkin fields with at most \(M\) spatial modes, max wave number \(K\),
viscosity \(\nu>0\), and energy \(E(t)\le E_0\),

\[
\Omega(t) \le \max\bigl(K^2 E_0,\; 6 M E_0^2 / \nu^2\bigr)
\qquad\forall t\ge 0.
\]

Uses \(\|\nabla\omega\|^2 \ge 2\Omega^2/E\) and the L-0002 embedding for \(\|\nabla u\|_\infty\).

## Caps (E0=0.5, ν=0.1)

| N | L-0001 (t=0.02) | L-0002 | L-0003 uniform |
|---|-----------------|--------|----------------|
| 8 | 8.779820e+01 | FAIL | 1.875000e+04 |
| 12 | 1.061457e+04 | FAIL | 5.145000e+04 |
| 16 | 1.203975e+11 | FAIL | 1.996500e+05 |

Still ≫ C-0001 numerical target (~11), but **closes globally in time** for each fixed N.
