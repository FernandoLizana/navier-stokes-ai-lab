# Function Spaces

**Purpose:** Fix the functional setting used in NS-MRL documentation and code comments.  
**Evidence:** N7 (standard definitions; no new theorems).

## 1. Periodic Sobolev spaces on \(\mathbb{T}^3\)

For \(s\in\mathbb{R}\),

\[
H^s(\mathbb{T}^3) = \Bigl\{ u=\sum_{k\in\mathbb{Z}^3}\hat u_k e^{ik\cdot x}:
\sum_{k}(1+|k|^2)^s|\hat u_k|^2 <\infty\Bigr\}.
\]

Divergence-free subspace:

\[
H^s_\sigma = \{ u\in H^s(\mathbb{T}^3)^3 : \nabla\cdot u = 0\}
\]
(mean-zero variants denoted with a dot when needed).

Smooth data for Route B live in \(\bigcap_s H^s_\sigma\).

## 2. Critical and subcritical spaces (literature pointers)

| Space | Role in NS theory |
|-------|-------------------|
| \(L^3\), \(\dot H^{1/2}\) | Energy-critical / scaling-critical |
| \(\mathrm{BMO}^{-1}\) | Koch–Tataru mild solutions (small data) |
| Besov \(\dot B^{-1+3/p}_{p,\infty}\) | Critical Besov frameworks |
| Gevrey / analytic classes | Analyticity radius; Foias–Temam type estimates |

NS-MRL exploratory diagnostics track Sobolev norms \(H^s\) for selected \(s\), enstrophy, and a crude analyticity-radius proxy from spectral decay. These diagnostics are **not** proofs of membership in a space for all time.

## 3. Weak / Leray–Hopf solutions

Leray–Hopf weak solutions satisfy the energy **inequality** and are known to exist globally for \(L^2\) data, but uniqueness and smoothness remain open in 3D. Partial regularity: Caffarelli–Kohn–Nirenberg (one-dimensional singular set in space-time for suitable weak solutions).

## 4. Computational truncation

Code works with Galerkin projections \(P_N\) onto modes \(|k|_\infty \le N/2\) (FFT grid). The continuum field is \(u = u_N + u_{\mathrm{tail}}\). Exploratory solvers ignore rigorous tail bounds (`≤ N2–N3`). Validated modules must bound \(u_{\mathrm{tail}}\) (`≥ N5` target).

## 5. Norm conventions in code

- \(L^2\) energy: \(E = \frac12 \sum_k |\hat u_k|^2\) with NS-MRL Fourier normalization (see `python/ns_exploration/spectral/fourier_conventions.py`).
- Enstrophy: \(\mathcal{E} = \frac12\|\omega\|_{L^2}^2\).
- All modules **must** import the single convention module; dual normalizations are forbidden.

## Primary references

- H. Triebel, *Theory of Function Spaces*.
- H. Bahouri, J.-Y. Chemin, R. Danchin, *Fourier Analysis and Nonlinear Partial Differential Equations*.
- H. Koch & D. Tataru, *Acta Math.* (2001).
- L. Caffarelli, R. Kohn, L. Nirenberg, *Comm. Pure Appl. Math.* (1982).
