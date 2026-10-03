# C-0008 Terminal-Weighted Audit (first delivery)

> N2 exploratory. **C-0008 remains open.** Not a proof.

## A. Identity / normalization

- direct cubic self-consistency: `0.000e+00`
- fullsym C stream vs dense: rel `0.000e+00`
- Stokes floor hi: `40.824623079` vs frozen `40.82462307930434`
- Band t=T regression vs L-0045: rel `0.000e+00`
- Full dealias sparse t=T: `25.925922025` vs L-0048 `25.925922025377965` (rel `0.000e+00`)
- Sparse vs dense {1..6}: `1.484e-16`

## B. L-0048 file map

See `reports/C0008_L0048_FILE_MAP.md`.

## Decision

Proceed to Phase D interval cover on band; then shell decomposition.
