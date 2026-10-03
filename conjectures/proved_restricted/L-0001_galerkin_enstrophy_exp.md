# L-0001 — Galerkin enstrophy exponential bound

**Status:** proved_restricted (finite Fourier–Galerkin only)  
**Evidence:** N7  
**Route tag:** B (setting), **not** a Clay contribution

## Statement

Under the hypotheses in `python/ns_exploration/conjectures/l0001_galerkin_bound.py`,

$$\Omega(t) \le K^2 E_0 \exp\bigl(2 K \sqrt{M} \sqrt{2 E_0}\, t\bigr).$$

## Constants on NS-MRL grids (dealiased)

| N | K | M | Ω(0.02) cap |
|---|-----|-----|-------------|
| 12 | 5.19615 | 343 | 1.061457e+04 |
| 16 | 8.66025 | 1331 | 1.203975e+11 |

Compare to numerical conjecture C-0001 with M≈10.9941 (explicit tight target, unproved).
L-0001 is rigorous but far looser; it shows an explicit finite-N ceiling exists constructively.

## Clay

None. K,M grow with resolution.
