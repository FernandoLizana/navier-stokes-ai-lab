# C-0008 — Phase D/E failure analysis (full-dealias)



**Date:** 2026-07-31  

**Phase D (Frobenius):** ~39 min sprint + manifest  

**Phase E (best):** **complete** — 87/87 shells via `upgrade_manifest_best_bounds.py` (~24h wall)



## Executive summary



| Quantity | Phase D (Frobenius) | **Phase E (best, measured)** |

|----------|---------------------|------------------------------|

| Best route | **D** | **D** |

| `I_term^hi` | **112.839** | **78.804** |

| `I_*^lo` | **1.2993** | 1.2993 |

| Gap vs `I_*` | ~87× | **~61×** |

| `Ω(T)^hi` (Route D) | **80.72** | **68.69** |

| `M_target` | 41.284 | 41.284 |

| Verifier | **FAIL** | **FAIL** |

| C-0008 status | **`exploring`** | **`exploring`** |



## Phase E — routes (full-dealias, N5, measured)



| Route | `I_term^hi` | `Ω(T)^hi` | Closes? |

|-------|-------------|-----------|---------|

| A shell integral | 170.22 | 101.01 | No |

| **D perturbation** | **78.80** | **68.69** | No (best) |

| B adaptive Lipschitz | 549.43 | 235.08 | No |



Improvement vs Frobenius Route D: **30.2%** on `I_term^hi` — **insufficient** for sharp M.



## Phase E — implementation notes



| Approach | Outcome |

|----------|---------|

| Combined pass 87 shells + 1-inf | OOM (`MemoryError` / worker kill) |

| **Per-shell 1-inf pass** | **Success** — checkpoint every 500k pairs |

| Band extrapolation (×0.516) | **Wrong** — full 1-inf often **worse** than Frobenius per shell |



Manifests:

- `shell_manifest_frobenius_full.json` (Phase D)

- **`shell_manifest_best_full.json`** (Phase E, authoritative for cert)



## Root cause



1. **Triangle + per-shell bounds** remain far from true `\|M(t)\|_{op}` at D=6748.

2. L-0048 float `C(T)≈25.93` is not a certified upper bound.

3. Even best Route D leaves `Ω(T)^hi` **> M_target** by ~27 units.



## Recommended next steps



- **Do not** promote C-0008 without verifier PASS.

- This certification **route is refuted** at N5 relaxation level.

- New structure needed: sparse/block op-norms, stretch (L-0027), or different temporal route.



## What was achieved



- First **reproducible full-dealias** N5 manifests (Frobenius + best).

- Rigorous Phase E completion with independent verifier (FAIL as expected).

- Honest negative: **C-0008 sharp target not closed** by Phase D/E pipeline.

