# Numerical Validation

Exploratory solvers (Python FFT, CUDA) produce **candidates** (`≤ N3`).

Validated computation must control: rounding, truncation, discretization, aliasing, Fourier tails, time error, residual, Lipschitz / inverse bounds, and produce machine-readable certificates.

**Sprint 1:** only a **finite-mode** interval module exists (Leray, energy identity on Galerkin truncation, short convolution). This does **not** validate the continuum PDE.
