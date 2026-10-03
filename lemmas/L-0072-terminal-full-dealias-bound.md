# L-0072 — Terminal full-dealias integral bound



**Status:** **route refuted** (Frobenius and best bounds) — C-0008 still **`exploring`**  

**Parent:** L-0071  

**Evidence:** N5 integral pipeline + verifier FAIL



## Result (2026-07-31)



Full-dealias manifests (D=6748, 87 shells), MPFR N5:



| Method | Route D `I_term^hi` | `I_*^lo ≈ 1.299` |

|--------|---------------------|------------------|

| Frobenius (Phase D) | 112.84 | FAIL |

| **min(F,1-inf) per shell (Phase E)** | **78.80** | **FAIL** |



Route A (best Phase E): 170.22. Route B: 549.43.



**Does not close** sharp target M ≈ 41.284 (`Ω(T)^hi` ≈ 68.69 > M).



Conclusion (spec §13): **this certification route cannot close C-0008**. Not a disproof of C-0008.



## Phase E (complete)



- 87 independent per-shell passes: `certified_L_op_one_inf_single_shell_hi`

- Checkpoint: `upgrade_best_checkpoint.json`, `upgrade_shell_ck/shell_NNN.json`

- Output: `experiments/terminal_weighted/shell_manifest_best_full.json`



On full dealias, 1-inf is **not uniformly tighter** than Frobenius per shell; `best` takes the minimum rigorously.



## Verifier



```bash

python verify_c0008_terminal_certificate.py \

  certificates/CERT-L0072-C0008-terminal-full-dealias.json \

  experiments/terminal_weighted/shell_manifest_best_full.json

```



**Expected:** FAIL (`I_term^hi >= I_*^lo`, `omega_T_hi >= M_target`).



C-0008 remains **`exploring`**.

