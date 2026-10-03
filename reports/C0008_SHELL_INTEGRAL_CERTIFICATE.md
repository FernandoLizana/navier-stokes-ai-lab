# C-0008 — Shell integral certificate (Route A)



## Formula (convention L-0071 / Phase D)



Per shell `r`, certified `B_r ≥ ||M_r||_op` (MPFR upward; **best** = min(Frobenius, 1-inf) in Phase E manifest).



\[

\int_0^T C_{\mathrm{term}}(t)\,dt

\le

\sum_{r \in \mathcal R}

B_r \,

\frac{1 - e^{-2\nu r T}}{2\nu r}.

\]



Weight factor `r` is included in `M_r` at `t=T`; temporal factor `e^{-2νr(T-t)}` integrated exactly.



## Threshold



\[

I_\star^{\mathrm{lo}} = 2\sqrt2\,(M_{\mathrm{target}} - \mathrm{StokesFloor}_{\mathrm{hi}})

\approx 1.2993127556607498.

\]



## Pilot band `{1..6}` (N5, not all-IC)



| Method | Route A `I_hi` | Closes `I_*`? |

|--------|----------------|---------------|

| Frobenius | 2.410 | Yes (band only) |

| **best** | **1.244** | Yes (band only) |



**Route D (band):** closes on band — **rejected** for C-0008 (`full_dealias=false`).



## Full-dealias (complete 2026-07-31)



| Method | Route A `I_hi` | Route D `I_hi` | Closes C-0008? |

|--------|----------------|----------------|----------------|

| Frobenius | 425.38 | 112.84 | No |

| **best (87/87)** | **170.22** | **78.80** | **No** |



Manifest: `experiments/terminal_weighted/shell_manifest_best_full.json`  

Certificate: `certificates/CERT-L0072-C0008-terminal-full-dealias.json`



Verifier recomputes Route A from manifest blocks (`recompute_integral.py`).



## Honesty



- Triangle inequality across shells.

- Per-shell bounds are not tight operator norms.

- Does **not** assume `C(t)` constant nor maximum at `t=T`.

