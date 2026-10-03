# C-0008 SOS / Triad Stretch Sprint

**Verdict:** `sos_subclass_only` — all-IC C <= C_dagger remains **open**.

## Thresholds

| Quantity | Value | vs C_dagger (9.5620) |
|----------|-------|------------------|
| C_emp (N2) | 0.009311 | yes (empirical) |
| C_fullsym all-IC | 25.9259 | **no** |
| C from triad R_star | 158.69 | **no** |
| Best SOS C_ub_hi | 0.8438224893075237 | subclass |

## L-0030 triad

- R_star = 3148
- Low-slab hybrid Omega(T) worst ~ 144.16
- min(hybrid, defect) worst ~ 51.00

R_★=3148, n_triads=4816686, M=3374. ‖N‖≤2√(R_★EΩ): at Ω=E=0.5 gives 56.11 vs Young 78.69. Cubic C≤2√2√R_★=158.7 (C_F=164.3, C_†=9.562). Low-slab hybrid Ω(T)≤144.2, min(hybrid,defect)≤51 (M=41.28); closes C-0007: False.

## SOS bands attempted

| Lemma | C_ub_hi | vs C_dagger |
|-------|---------|-------------|
| CERT-L0052-sos-gram-band123-N24 | 0.8438224893075237 | yes |
| CERT-L0053-sos-gram-onepol-cr0008-N24 | 3.3560651946371634 | yes |
| CERT-L0054-sos-gram-band1234-N24 | 0.8940401324113961 | no |
| CERT-L0058-sos-gram-band12345-N24 | 1.3052615184195027 | yes |

## Next targets

1. Weighted triad / SOS flattening (L-0030 notes)
2. Extend N5 SOS one-pol greedy beyond C-R-0014
3. Hybrid min(triad, defect, SOS) on low slab ODE

Low slab Ω0≤Ω★=18.874042: C_†=9.56201 ⇒ Ω(T)≤41.284 (M=41.284). C_emp≈0.00931135; emp_closes=True. All-IC cubic at C=0 already has floor=46.2846>M (high slab needs L-0024/L-0026). Proved C≤C_†: False.
