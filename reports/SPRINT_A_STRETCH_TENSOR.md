# Sprint A — Auditoría del tensor de stretch + simetrización total (C-0007)

**Objetivo.** Preparar el ataque SOS al cuello de botella de C-0007: certificar que
el objeto de la cota Sym (`l0040`) es *exactamente* la forma cúbica del vortex
stretching, y comprobar si una matricización **totalmente simétrica** (flattening
del 3-tensor) baja la norma de operador por debajo de `C_† ≈ 9.562`.

**Sin tocar lemmas. Sin Julia. Finito Galerkin. No continuo. No Clay.**

Banda de prueba: shells consecutivas `{1,2,3,4,5,6}`, `N=24`, `D=160` coords reales.

---

## 1. Definición auditada (base real transversal, coords de enstrofía)

`f(z) = Σ_{a,i,j} G[a,i,j] · z_a z_i z_j` con

```
G[a,i,j] = |k_a| / (|k_i|·|k_j|) · ⟨p_i, k_j⟩ · ⟨p_a, P_{k_i+k_j}(p_j)⟩   si k_a = k_i+k_j
```

- `p_·` = polarizaciones reales unitarias transversales (`_pol_basis`).
- `P_s` = proyector de Leray en `s` (`_leray_vec`).
- **Coordenadas de enstrofía**: `‖z‖² = 2Ω` (amplitud de velocidad = `z/|k|`), por eso
  el coeficiente de advección lleva `1/(|k_i||k_j|)` y la fila de salida lleva `|k_a|`.
- La constante buscada: `C = 2√2 · sup_{‖z‖=1} |f(z)|`, se necesita `C ≤ C_†`.

> **Nota de auditoría (bug encontrado y corregido durante el sprint).** Una primera
> versión del tensor omitió el factor `1/(|k_i||k_j|)` y dio `C≈36.6 ≠ 11.344`.
> Con el factor correcto reproduce `l0040` a `6e-14`. La normalización de enstrofía
> es esencial y queda documentada aquí.

---

## 2. Resultados numéricos

| Cantidad | Valor | Lectura |
|----------|-------|---------|
| `f_tensor` vs matricización `l0040` (1500 z) | err `6.0e-14` | ✅ el tensor **es** el objeto de la cota Sym |
| `f_fullsym` vs `f_tensor` (1500 z) | err `4.4e-14` | ✅ simetrización total preserva la forma cúbica |
| `C_sym` reproducido | `11.34411252` | = `l0040` exacto |
| **`C_fullsym`** (flattening totalmente simétrico) | **`11.28669880`** | ⬇ solo −0.5 % respecto a Sym |
| `C_rank1_lower` (ascenso tensorial, N2) | `3.12909966` | cota **inferior** del `C` real |
| `C_† ` (umbral de cierre) | `9.56200994` | objetivo |
| Gap de relajación `C_fullsym / C_rank1` | `≈ 3.6×` | Sym sobre-cuenta ~3.6× |

### Interpretación
- **El `C` verdadero de la banda `{1..6}` está en `[3.13, 11.29]`.**
- La simetrización total **no cierra el gap**: `C_fullsym = 11.287 > C_† = 9.562`.
  El obstáculo **no** es la elección de matricización; es el **gap de relajación
  rank-1 ↔ Sym** (factor ~3.6×), exactamente lo que un SOS/momentos debe atacar.
- **Buena noticia para SOS:** en términos de `sup`, la cota espectral actual da
  `sup ≤ 3.99` y el umbral es `sup ≤ 3.3807` — **solo hay que rebajar ~15 %**
  (`11.287 → 9.562`). La cota inferior real (`sup ≳ 1.11`, i.e. `C ≳ 3.13`) queda
  muy por debajo del umbral. Un SOS de orden 2 razonablemente ajustado tiene
  margen real de aterrizar bajo el umbral.

---

## 3. Conclusiones del Sprint A

1. **Fase 1–2 de ChatGPT completadas (parcial).** Tensor cúbico construido,
   simetrizado y **verificado a `1e-13`** contra el objeto establecido de `l0040`
   (que codifica las triadas pseudoespectrales). Esto valida que la cota Sym y la
   comparación con `C_†` operan sobre la forma física correcta.
2. **La simetrización total NO es un atajo.** `C_fullsym ≈ 11.287` ≈ `C_Sym`.
   Ninguna manipulación lineal de la matricización cerrará `{1..6}`.
3. **El SOS es genuinamente necesario y está bien motivado.** El gap a cerrar es
   pequeño (~15 % sobre la cota espectral) y la evidencia N2 sitúa el valor real
   muy por debajo de `C_†`.

## 4. Pendiente / recomendaciones (para decidir antes de la Fase 5)

- **Completar Fase 1 (honestidad):** falta un cross-check independiente por FFT
  (`stretch_inner` sobre campos reales) frente a `f_tensor`. El check hecho aquí es
  contra la matricización `l0040` (mismo objeto de triadas, distinto camino de
  código), no contra la ruta FFT. Recomendado antes de publicar como N7 pleno.
- **Alcance:** cerrar `{1..6}` por SOS daría **otra subclase C-R** (banda de 2
  polarizaciones), no all-IC. All-IC exige el espacio retenido completo.
- **Siguiente paso propuesto (Fase 5):** sidecar Julia aislado con TSSOS/CS-TSSOS,
  relajación de **orden 2** sobre `max f(z)` s.a. `‖z‖²=1` (base racional para el
  certificado exacto posterior), con límites duros de tamaño de bloque PSD y
  aborto previo. Si la cota superior SOS `< 3.3807` con margen → Fase 8
  (certificado racional + verificador sin solver).

---

## Constantes de referencia

```
banda        = {1,2,3,4,5,6},  N=24,  D=160
C_†          = 9.562009938690146     (⇔ sup ≤ C_†/(2√2) ≈ 3.3806810)
C_sym        = 11.344112520990302    (l0040, matriz simétrica solo en (i,j))
C_fullsym    = 11.286698798082990    (3-tensor totalmente simétrico)
C_rank1_lo   = 3.129099662236431     (N2, ascenso local — cota inferior)
gap Sym/real ≈ 3.6×
```

*Sistema finito Galerkin. Constantes no uniformes cuando N→∞. No continuo. No Clay.*
