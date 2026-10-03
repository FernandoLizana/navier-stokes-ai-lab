# Sprint — L-0044 physical Hermitian fullsym + C-R-0009

| Item | Value |
|------|-------|
| **C-R-0009** (Hermitian shells {1..6}) | C_fullsym≈2.993326 ≤ C_† |
| (i,j)-only Sym (same tensor) | ≈11.344113 > C_† |
| FFT ↔ tensor max rel err | 5.15e-14 |
| C-0007 all-IC | still open |

Key audit: physical triad uses factor `i`; full symmetrization then closes {1..6}.

FINITE Galerkin only. Not continuum. Not Clay.
