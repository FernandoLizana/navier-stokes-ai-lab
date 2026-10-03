# L-0027 Stretch Sprint — Post-Terminal Repair (2026-08-06)

**C-0008:** `exploring`

## Results

| Check | Value | Closes? |
|-------|-------|---------|
| C_dagger (techo) | ~9.56 | — |
| C_emp_max (low-slab random) | ~0.0093 | low-slab only |
| C_fullsym all-IC | ~25.9 | **NO** vs C_dagger |
| Greedy SOS extrap | ~5.39 | subclass |
| Terminal route (legacy) | **REFUTED** | repair supersedes |

## Viable paths (sprint output)

1. **L-0027 stretch** — prove all-IC `C ≤ C_dagger` (fullsym ~25.9 refutes naive Shor)
2. N5 SOS subclasses C-R-0002..0014 (proved restricted)
3. Weighted triad / SOS flattening
4. Non-partition terminal (L-0074 refuted finer clusters)

## Relation to repair baseline

Repair Route A: I_hi ≈ **170.59** (honest, verify PASS).  
Repair cluster L-0073R: I_hi ≈ **38.52** (audit PASS, sketch).  
L-0027 does **not** bypass terminal bounds — structural all-IC path remains orthogonal.

Legacy L-0073 ~37.92 and terminal Phase D/E certs remain **refuted** by repair verifiers.
