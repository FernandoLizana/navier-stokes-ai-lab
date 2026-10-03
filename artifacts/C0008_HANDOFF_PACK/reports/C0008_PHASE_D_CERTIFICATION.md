# C-0008 — Phase D/E certification report (complete)



**Date:** 2026-07-31  

**Status:** C-0008 **`exploring`** (verifier FAIL — honest negative result)



---



## A. ¿Qué cantidad se certificó?



\[

I_{\mathrm{term}} := \int_0^T C_{\mathrm{term}}(t)\,dt

\]



con `C_term(t) = 2√2 ‖M(t)‖_{op}` y pesos terminales L-0071 en fullsym Galerkin.



## B. ¿Full-dealias?



**Sí** — 87 shells dealias, D=6748.



| Manifest | Método | Archivo |

|----------|--------|---------|

| Frobenius | Phase D sprint | `shell_manifest_frobenius_full.json` |

| **Best** | Phase E (87/87) | **`shell_manifest_best_full.json`** |



## C. ¿N cubiertos?



**N7 sketch:** wavevectors n=12,16,20 ⊂ n=24. Certificado en N=24 bajo monotonicidad de modo (no N6 formal).



## D. Backend riguroso



**gmpy2 MPFR 2.3.1**, RoundUp/RoundDown, 128 bits.



## E. Acotación de norma



- Frobenius MPFR upward (N5)

- `√(‖M‖₁‖M‖∞)` MPFR upward (N5), **por shell** (pass D²/memoria acotada)

- `min(Frobenius, 1-inf)` por shell en manifest **best**



## F. `I_term^hi` (full-dealias)



### Phase D — Frobenius



| Route | `I_term^hi` | `Ω(T)^hi` |

|-------|-------------|-----------|

| A | 425.38 | 191.22 |

| **D** | **112.84** | **80.72** |

| B | 1378.55 | 528.22 |



### Phase E — best (87/87 shells, 2026-07-31)



| Route | `I_term^hi` | `Ω(T)^hi` | Closes C-0008? |

|-------|-------------|-----------|----------------|

| A | **170.22** | 101.01 | No |

| **D** | **78.80** | **68.69** | **No** |

| B | 549.43 | 235.08 | No |



Mejora Route D vs Frobenius: **~30%** en `I_term^hi`; sigue **~61×** por encima de `I_*^lo`.



Band pilot (D=160): Route A best **1.244** — cierra `I_*` en banda solo (`full_dealias=false`).



## G. `I_*^lo`



**1.2993127556607498** (MPFR recalculado)



## H. Margen (Route D best)



**Negativo** — `I_term^hi - I_*^lo ≈ +77.5`.



## I. `Ω(T)^hi`



**68.69** (Route D best) > M_target **41.284**.



## J. Verificador



**FAIL** (esperado):



```powershell

python verify_c0008_terminal_certificate.py `

  certificates/CERT-L0072-C0008-terminal-full-dealias.json `

  experiments/terminal_weighted/shell_manifest_best_full.json

```



## K. Niveles de evidencia



| Pieza | Nivel |

|-------|-------|

| L-0071 identidad | N6 sketch |

| L-0048 float | N2 (no upper) |

| Frobenius / 1-inf per shell | N5 |

| N≤24 | N7 |

| C-0008 sharp | **abierto** |



## L. Limitaciones



1. Triangular + per-shell bounds siguen muy laxos en D=6748.

2. En full-dealias, 1-inf **no siempre** baja Frobenius por shell.

3. Route D no demuestra monotonicidad en t.

4. **Esta ruta no cierra C-0008** (spec §13: refutación de ruta, no de C-0008).



---



## Reproducir



```powershell

# Verificar artefactos precomputados (~segundos)

.\reproduce_c0008_certificate.ps1



# Recomputar todo (~24h+ Phase E)

.\reproduce_c0008_certificate.ps1 -Full

```



Artefactos: `certificates/CERT-L0072-C0008-terminal-full-dealias.json`, `experiments/terminal_weighted/phase_d_run.json`, `lemmas/L-0072-terminal-full-dealias-bound.md`.

