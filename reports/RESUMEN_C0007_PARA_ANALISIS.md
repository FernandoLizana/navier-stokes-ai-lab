# NS-MRL Ruta B — Resumen técnico para análisis externo (C-0007)

> **Propósito de este documento.** Es un resumen **autocontenido** del trabajo hecho
> en un proyecto de cotas de enstrofía para Navier–Stokes **truncado (Galerkin/pseudoespectral)**.
> Está pensado para que un revisor externo (p. ej. ChatGPT) lo analice y opine sobre
> viabilidad y siguientes pasos. **No requiere acceso al repositorio.**

---

## 0. Disciplina de alcance (leer primero)

- **Todo es finito-dimensional.** Cada enunciado es sobre una truncación Galerkin /
  pseudoespectral en una malla `N³` con dealias 2/3. **Nada de esto es una afirmación
  sobre el PDE continuo de Navier–Stokes, y nada tiene implicación con el problema Clay.**
- Las constantes de cierre escalan como `∼ M/ν²` y **divergen cuando `N→∞`** (`M→∞`),
  por lo que **no** hay puente al continuo por esta vía. Es un ejercicio de
  demostración→certificado→verificador sobre enunciados finitos que sí se pueden cerrar.
- Escala de evidencia: **N2** = simulación en float; **N5** = constante certificada por
  intervalos; **N6** = conjetura; **N7** = demostración a mano del sistema truncado.

---

## 1. Setup y parámetros

- Dominio `T³ = [0,2π)³`, campo incompresible `u` (`div u = 0`).
- Energía `E(t) = ½‖u‖₂²`, enstrofía `Ω(t) = ½‖∇u‖₂² = ½‖ω‖₂²` (`ω = curl u`).
- Parámetros fijos del objetivo: **`ν = 0.1`, `E₀ = 0.5`, `T = 0.02`, `N ≤ 24`, dealias 2/3.**
- Energía no creciente (`ν ≥ 0`) ⇒ `E(t) ≤ E₀`.
- Base: modos de Fourier retenidos tras dealias, con proyección de Leray
  (2 polarizaciones transversales por vector de onda `k`). `|k|² ∈ ℤ` indexa "shells".

---

## 2. El objetivo abierto: C-0007

> **C-0007 (all-IC, N6, abierta).** Para *todo* campo incompresible de Fourier en `N ≤ 24`
> con `E ≤ 0.5`, la NS Galerkin dealias (`ν = 0.1`) cumple
> **`Ω(0.02) ≤ M := 41.28399951`.**

El objetivo vive en un **hueco numérico estrecho**:

| Cota | Valor | Rol |
|------|-------|-----|
| Stokes floor | `40.82462308` | máximo del sistema **lineal** de Stokes en `T` (no acota la NS no lineal) |
| **M (objetivo C-0007)** | **`41.28399951`** | `≈ 1.011 × floor` |
| Majorante L-0024 | `41.74337594` | mejor cota all-IC probada hoy (cierra C-0005/0006, **no** C-0007) |
| Envelope cap | `73.5` | cota cruda superior |

Es decir: L-0024 ya prueba `Ω(0.02) ≤ 41.743` para todo IC, pero **no baja hasta 41.284**.

---

## 3. Estrategia de dos franjas ("slabs") sobre Ω₀

Se parte el espacio de datos iniciales según la enstrofía inicial `Ω₀`, con umbral
`Ω★ ≈ 18.874`:

### 3.1 Franja alta `Ω₀ ≥ Ω★` — **CERRADA**
El ODE de defecto espectral (L-0024) ya lleva `Ω(T) ≤ M`. Empaquetado como **C-R-0002**
(vía L-0024 + L-0026). Sin pendientes.

### 3.2 Franja baja `Ω₀ ≤ Ω★` — **ABIERTA (cuello de botella)**
Aquí el ODE cúbico (L-0023/L-0027) cierra `Ω(T) ≤ M` **si y solo si** el término de
*stretch* (vortex stretching) satisface una cota cúbica

```
|⟨ω, curl N(u)⟩| ≤ C · Ω^{3/2}     con     C ≤ C_† ≈ 9.562009939.
```

- Numéricamente (N2), el `C` real observado es **`≪ C_†`** (órdenes de 0.1–1.3).
- Pero **`C ≤ C_†` no está probado all-IC.** Ese es el único hueco lógico que impide C-0007.

---

## 4. El problema matemático exacto del cuello de botella

Tras proyección de Leray y matricización, el stretch es una **forma cúbica** en los
coeficientes `z ∈ ℝ^D` del campo (dimensión `D` = nº de modos×polarizaciones retenidos):

```
stretch(z) = z · L(z zᵀ)
```

donde `L` es un operador lineal (tensor de triadas) que se puede montar exactamente.
La constante buscada es

```
C = 2√2 · sup_{‖z‖=1} | z · L(z zᵀ) |,      hace falta   C ≤ C_†  ⇔  sup ≤ C_†/(2√2) ≈ 3.3807.
```

### 4.1 La relajación "Sym" y el GAP
El único majorante **probado** que tenemos reemplaza el tensor rank-1 `zzᵀ` por una matriz
simétrica **general** `Z` con `‖Z‖_F = 1`:

```
|stretch| ≤ ‖L|_Sym‖_op · ‖z‖³     ⇒     C_Sym = 2√2 · ‖L|_Sym‖_op.
```

Esto es un **SVD denso exacto (N7)**, pero **sobre-cuenta** porque admite `Z` que no son
`zzᵀ` (rank > 1). El gap resultante es enorme:

| Objeto | Valor de C | ¿≤ C_†? |
|--------|-----------|---------|
| `C_Sym` en shells `{1..6}` (2 polarizaciones) | `≈ 11.344` | ❌ (falla por ~19%) |
| `C` real rank-1 (empírico, mismas shells) | `≈ 0.23 – 1.3` | ✅ (con enorme margen) |

**Conclusión clave:** la barrera **no es física** (el stretch real es pequeñísimo);
es el **gap de relajación** de Sym. Certificar el máximo **rank-1** (cuártica no convexa)
es exactamente un problema **SOS / de momentos (Positivstellensatz)**.

### 4.2 Todo lo que se puede hacer SIN un solver SOS ya se agotó
Con solo `numpy`/`scipy` (sin `cvxpy`/`MOSEK`/`SCS`) se probaron majorantes alternativos y
**todos** quedan por encima de `C_†` en la banda consecutiva `{1..6}`:

| Técnica (techo) | C obtenido | Lemma |
|-----------------|-----------|-------|
| Sym relajación (SVD denso) | `≈ 11.344` | L-0040/41 |
| Column-triangle majorant | `≈ 80` | L-0042 |
| A/B block hybrid (`A={1..5}, B={6}`) | `≈ 16.21` | L-0042 |
| Shell-pair→target block (muestreado, cota inf.) | `≳ 14.14` | L-0043 |
| Ky-Fan / suma SVD (4 términos) | `≳ 12.65` | L-0043 |
| Rank-1 Rayleigh (cota **inferior**, no certifica) | `≈ 0.23` | L-0042 |

---

## 5. Resultados POSITIVOS: subclases de C-0007 ya probadas (N7)

Cada `C-R-XXXX` es un subconjunto de datos iniciales donde **sí** se prueba `Ω(0.02) ≤ M`.
Todas usan `L-0026` (franja alta) + `L-0027` (franja baja con `C ≤ C_†`).

| ID | Subclase (dominio de IC) | Constante | Lemma |
|----|--------------------------|-----------|-------|
| **C-R-0001** | IC mono-radial (una shell) | M≈71.73 | L-0020+L-0021 |
| **C-R-0002** | `Ω₀ ≥ Ω★ ≈ 18.87` (franja alta, all-IC) | — | L-0024+L-0026 |
| **C-R-0003** | mono-radial con `α_r = 0` (incl. shell 147) | — | Stokes–Duhamel |
| **C-R-0004** | soporte de Fourier en shells `{1,2,3}` | — | L-0038 |
| **C-R-0005** | soporte en 14 shells (incl. `{1,2,3,147,…}`) | C≈9.519 | L-0039 |
| **C-R-0006** | soporte en shells `{1,2,3,4,5}` | C_Sym≈9.015 | L-0040 |
| **C-R-0007** | soporte en 8 shells con r=6: `{1,2,3,4,6,11,12,25}` | C_Sym≈9.414 | L-0041 |
| **C-R-0008** | **una sola polarización** por modo, shells `{1,2,3,4,5,6,8}` | C≈8.476 | L-0043 |

**Hallazgo reciente (C-R-0008):** al congelar **una polarización** por vector de onda
(en vez de las 2), Sym cierra una banda consecutiva más larga (`{1..8}`, `C≈8.48 ≤ C_†`);
añadir la shell 9 vuelve a fallar (`C≈9.85 > C_†`).

**Limitación estructural:** estas subclases son familias muy particulares (soportes
espectrales fijos o 1 polarización). **No se acumulan a all-IC**: enumerarlas nunca cierra
C-0007 por sí solo.

---

## 6. Inventario de lemmas puente

- **L-0001…L-0006** — cotas cúbicas base del sistema truncado (crudo, viscoso, uniforme-en-tiempo, refinamiento por shells).
- **L-0007…L-0017** — puente shell-IC → malla completa; sub-conjetura C-S-0002 a `T=0.02` **sigue abierta** en horizonte largo.
- **L-0018** — envelope; cierra C-0003/C-0004/C-S-0002 (envelope 37.5).
- **L-0019** — Stokes floor. **L-0020/L-0021** — Duhamel H¹ + cuártica radial (C-R-0001).
- **L-0022/L-0023** — bootstrap γ + disipación cúbica (el ODE de la franja baja).
- **L-0024** — defecto espectral ⇒ **C-0005** (`M≈61.24`); majorante all-IC `≈41.743`.
- **L-0025** — defecto multi-N ⇒ **C-0006** (N=32, `M≈67.77`).
- **L-0026** — franja alta `Ω₀≥Ω★` cumple M (⇒ C-R-0002); franja baja es techo estructural.
- **L-0027** — ODE cúbico condicional en franja baja: cierra si `C ≤ C_†≈9.562`.
- **L-0028…L-0037** — **techos** (embeddings, Frobenius, triadas `R★`, `S_max`, Shor firmado,
  diferencias de shell): todos dan `C ≫ C_†`; ninguno cierra la franja baja all-IC.
- **L-0038/L-0039** — Shor por SVD de bandas / soportes dispersos ⇒ C-R-0004/0005.
- **L-0040/L-0041** — Shor Sym-restringido ⇒ C-R-0006/0007; techo consecutivo `{1..6}`.
- **L-0042** — intentos de superar Sym en `{1..6}` (triangle/hybrid/rank-1): todos fallan.
- **L-0043** — one-pol Sym ⇒ **C-R-0008** (`{1..8}`); techos shell-block/Ky-Fan.

### Conjeturas rechazadas (por completitud)
- **C-0001, C-0002** — refutadas (modo de Stokes `|k|²=75` da `Ω(0.02)≈27.78 > M≈20.32`).
- **C-S-0001** — refutada (N2).

---

## 7. Estado y pendientes

### Cerrado
- **C-0003, C-0004, C-0005, C-0006** (all-IC, distintos M) — probadas.
- **C-R-0001 … C-R-0008** — subclases de C-0007 probadas.
- Franja alta de C-0007 (`Ω₀ ≥ Ω★`) — cerrada (C-R-0002).

### Abierto (bloqueante)
1. **C-0007 all-IC** — falta la franja baja `Ω₀ ≤ Ω★`.
2. El único hueco lógico: probar `C ≤ C_† ≈ 9.562` para el stretch **rank-1** (`z zᵀ`),
   equivalente a certificar `sup_{‖z‖=1} |z·L(z zᵀ)| ≤ 3.3807`.
3. Es un **problema SOS / de momentos de grado 4, no convexo**. Todo majorante convexo
   "barato" (Sym, triangle, hybrid, shell-block, Ky-Fan) sobre-cuenta y supera `C_†`.
4. Herramienta faltante: un **solver SOS** (no hay `cvxpy`/`MOSEK`/`SCS` instalados).

### Diferido (infra)
- Lean N8 (esqueleto hecho; ODE/PDE pendiente), CAP/Newton–Kantorovich, Julia validated
  numerics, publicación parcial. Evaluación Clay: **solo tras una técnica continua general**
  (no aplicable a esta vía finita).

---

## 8. Preguntas concretas para el revisor externo

1. **¿El gap de relajación es cerrable con SOS de grado bajo?** El objeto es
   `p(z) = C_†²·‖z‖⁶ − (2√2)²·(z·L(z zᵀ))²`; ¿es plausible una descomposición
   `p = σ₀ + Σ σ_i·g_i` (con `g_i` las restricciones de div-free/energía) de grado 4–6
   que sea numéricamente tratable en `D ~ 80–160`?
2. **¿Hay una cota rank-1 exacta** que evite SOS? P. ej. explotar la estructura de triada
   (`k = p+q`), antisimetría del stretch, o identidades de conservación (`∑ energía`,
   helicidad) para acotar `z·L(z zᵀ)` directamente por debajo de `3.3807`.
3. **¿Merece la pena?** Cerrar C-0007 es un teorema finito estrecho (`Ω(0.02)≤41.284` en
   `N≤24`), **sin** valor Clay. ¿Es mejor (a) invertir en el solver SOS, (b) seguir el
   juego infinito de subclases C-R, o (c) declarar el estado actual como conclusión honesta?
4. **Sanidad del target M.** `M=41.284` está entre floor `40.825` y majorante L-0024
   `41.743`. ¿El valor exacto de `M` (`≈1.011×floor`) es defendible, o conviene relajarlo
   hacia `41.743` para tener un enunciado all-IC ya demostrado (= L-0024) y cerrar de una vez
   con un M ligeramente mayor?

---

## 9. Constantes numéricas de referencia (para verificación)

```
ν            = 0.1
E₀           = 0.5
T            = 0.02
N_max        = 24        (dealias 2/3)
M (C-0007)   = 41.283999509509286
Stokes floor = 40.82462307930434
L-0024 maj.  = 41.74337593971423
envelope cap = 73.5
Ω★ (umbral)  = 18.874042281544565
C_†          = 9.562009938690146      (⇔ sup rank-1 ≤ C_†/(2√2) ≈ 3.3806810)
C_Sym {1..6} = 11.344112520990302     (2 polarizaciones, > C_†)
C   {1..8}   = 8.476437931112338      (one-pol, C-R-0008, ≤ C_†)
```

*Sistema finito Galerkin. No continuo. No Clay.*
