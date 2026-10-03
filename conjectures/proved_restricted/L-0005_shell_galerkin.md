# L-0005 — Shell-truncated Galerkin enstrophy bound

**Status:** proved_restricted  
**Evidence:** N7 (via L-0001 + L-0004)  
**Clay:** none

## Statement

Let \(P_{K_0}\) be the Fourier projector onto modes with \(|k|\le K_0\) (and dealias mask).
For the Galerkin ODE \(\partial_t u = P_{K_0}[-P(u\cdot\nabla)u + \nu\Delta u]\),
the support remains in the shell and

\[
\Omega(t)\le \min\bigl(\text{L-0001}(K_0,M_{K_0},t),\;\text{L-0003}(K_0,M_{K_0})\bigr).
\]

Example (E0=0.5, ν=0.1, t=0.02, N=12): \(K_0=4\) gives a finite explicit ceiling
(see `L-0005.json`).

## Relation to C-S-0001

C-S-0001 allows **full** dealias evolution from shell ICs (higher modes may appear).
L-0005 does **not** prove C-S-0001; it proves the truncated companion problem.
