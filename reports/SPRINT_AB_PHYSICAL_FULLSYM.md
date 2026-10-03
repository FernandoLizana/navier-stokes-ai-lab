# Sprint A+B — Tensor audit, FFT cross-check, L-0044 / C-R-0009

**Finito Galerkin. No continuo. No Clay.**

## Hallazgo central

El tensor usado en L-0040 (amplitudes reales por modo complejo, **sin** factor `i`)
**no** es el stretch físico `⟨ω, curl N⟩`. El objeto físico correcto es:

```
N_s = -P_s Σ_{p+q=s} i (û_p · q) û_q
stretch = ⟨ω, curl N⟩
```

en la base Hermitiana real `(c,s)` por modo×pol del semi-espacio, con `‖z‖² = 2Ω`.

Sobre ese tensor, la **simetrización total** cambia todo:

| Objeto | C | ¿≤ C_†≈9.562? |
|--------|---|----------------|
| Sym `(i,j)` físico `{1..6}` | 11.344 | ❌ |
| **Fullsym físico `{1..6}`** | **2.993** | ✅ |
| Rank-1 lo (N2) | 0.542 | (cota inferior) |
| FFT ↔ tensor (200 z) | err **2.7e-14** | ✅ |

**C-R-0009** = campos reales con soporte Fourier en shells `{1,2,3,4,5,6}` ⇒ `Ω(0.02)≤M`.

## Exploratorio (mismo majorante)

| Bandas | D | C_fullsym | ok |
|--------|---|-----------|-----|
| r≤10 | 292 | 4.26 | ✅ |
| r≤14 | 500 | 6.14 | ✅ |
| r≤20 | 776 | 7.70 | ✅ |
| r≤24+ | — | OOM (G denso) | — |

All-IC sigue abierto: hace falta full-mask sin almacenar `G` denso (o SOS sparse).

## Archivos

- `experiments/sprintA_stretch_tensor.py` — tensor Galerkin / fullsym A
- `experiments/sprintB_physical_tensor.py` — tensor Hermitiano + FFT
- `conjectures/l0044_physical_fullsym.py` — L-0044 + C-R-0009
- `tools/sos_julia/` — scaffold (Julia no instalado aún)

## Errores de escala documentados

1. Omitir `1/(|k_i||k_j|)` en coords de enstrofía → C≈36 ≠ 11.34.
2. Omitir factor `i` / no imponer Hermiticidad → stretch FFT ≈ 0 o desacuerdo total.
3. Sym solo en `(i,j)` sobre el tensor físico **no** cierra; hace falta fullsym.
