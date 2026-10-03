# NS-MRL Ruta B — Paquete completo para revisión externa (ChatGPT / revisor)

> **Actualizado:** post L-0070 (Opción A). **Autocontenido:** no requiere repo.
> **Alcance:** Galerkin/pseudoespectral finito T³, N≤24, dealias 2/3. **NO continuo. NO Clay.**

---

## A. Qué pedirle al revisor

1. ¿Hay una **nueva dirección matemática** para cerrar **C-0008** (M≈41.284) sin repetir la escalera Sym/SOS que ya falló en all-IC?
2. ¿El gap **rank-1 vs Sym** es cerrable con otra técnica (no SOS brute)?
3. ¿Merece la pena seguir atacando C-0008, o el estado actual (C-0007 all-IC a 41.743) es conclusión honesta?
4. ¿Qué estructura espectral / física (triadas, helicidad, cancelaciones) podría bajar el majorante L-0024 ~0.46 unidades?

---

## B. Setup fijo

| Parámetro | Valor |
|-----------|-------|
| Dominio | T³ = [0,2π)³, u incompresible |
| ν | 0.1 |
| E₀ (energía máx IC) | 0.5 |
| T | 0.02 |
| N máx | 24 (dealias 2/3) |
| K² (shell máx \|k\|²) | 147 (8 modos ±7,±7,±7) |
| gap (K² − r_next) | 13 |
| m_K (multiplicidad shell 147) | 8 |
| envelope cap K²E₀ | 73.5 |

**Energía** E = ½‖u‖², **enstrofía** Ω = ½‖∇u‖² = ½‖ω‖².

---

## C. Estado actual de conjeturas (2026-07-28)

### CERRADO (all-IC, N7)

| ID | Enunciado | M / cota | Método |
|----|-----------|----------|--------|
| **C-0003** | all-IC envelope | 37.5 | L-0018 |
| **C-0004** | all-IC | ~40.8 | L-0018 |
| **C-0005** | all-IC | M≈61.24 (majorante real ≈41.743) | L-0024 |
| **C-0006** | all-IC N≤32 | M≈67.77 | L-0025 |
| **C-0007** | **all-IC** | **Ω(0.02) ≤ 41.74337594** | **L-0070 / L-0024** |

### ABIERTO

| ID | Enunciado | Target M | Bloqueo |
|----|-----------|----------|---------|
| **C-0008** | sharp refinement de C-0007 | **41.28399951** | Gap 0.459 vs L-0024; necesita C≤C_† en low-slab o argumento nuevo |

### RECHAZADO (refutado)

- **C-0001, C-0002** — modo Stokes |k|²=75 da Ω(T)≈27.78 > M propuesto
- **C-S-0001** — refutado N2

---

## D. Constantes numéricas clave

```
Stokes floor          = 40.82462307930434
C-0008 sharp M        = 41.283999509509286   (= mid(floor, L-0024 maj))
C-0007 proved M       = 41.74337593971423    (= L-0024 majorant)
Gap sharp vs proved   = 0.459376430204944
L-0024 worst Ω0       = 0.5 (= E0, Stokes-like)
Ω★ (umbral slab)      = 18.874042281544565
C_† (stretch target)  = 9.562009938690146
  ⇔ sup rank-1        ≤ C_†/(2√2) ≈ 3.3806810
C_emp (N2, max seen)  ≈ 0.0437  (≪ C_†, no prueba)
D full dealias        = 6748 modos (one-pol ~ half)
```

---

## E. Estrategia de dos franjas (L-0026 / L-0027)

Partir IC por enstrofía inicial Ω₀ con umbral Ω★≈18.874:

### E.1 Franja alta Ω₀ ≥ Ω★ — CERRADA (C-R-0002)
ODE defecto espectral L-0024 ⇒ Ω(T) ≤ M sharp (41.284).

### E.2 Franja baja Ω₀ ∈ [E₀, Ω★] — ABIERTA para C-0008
ODE cúbico L-0023/L-0027 cierra Ω(T)≤M **si y solo si**

```
|⟨ω, curl N(u)⟩| ≤ C · Ω^{3/2}    con    C ≤ C_† ≈ 9.562
```

- **C ≤ C_† NO está probado all-IC.**
- C_emp ≈ 0.04 (simulación); margen enorme pero no certificado.
- ODE cúbico solo (C=0) da floor≈46.28 > M → **no cierra all-IC sin L-0024**.

### E.3 C-0007 all-IC (Opción A, L-0070)
**Sin case-split:** L-0024 majorant uniforme, Ω(T) peor en Ω₀=E₀, monótona no creciente en Ω₀.
⇒ **Ω(0.02) ≤ 41.743 para todo IC.** No usa stretch C_†.

---

## F. Problema matemático del cuello de botella (stretch rank-1)

Coeficientes z ∈ ℝ^D (modos dealias × polarizaciones Leray):

```
stretch(z) = z · L(z zᵀ)
C = 2√2 · sup_{‖z‖=1} |z · L(z zᵀ)|
```

Necesario para C-0008 vía L-0027: **C ≤ 9.562**.

### Relajación Sym (SVD denso, N7) — GAP ENORME

```
|stretch| ≤ ‖L|_Sym‖_op · ‖z‖³   ⇒   C_Sym = 2√2 · ‖L|_Sym‖_op
```

Sym permite Z simétrica general, no solo zzᵀ (rank 1). **Sobre-cuenta.**

| Soporte / método | C_Shor o C_Sym | ¿≤ C_†? |
|------------------|----------------|---------|
| shells {1..6}, 2 pol | ≈11.34 | ❌ |
| one-pol {1..8} C-R-0008 | 8.476 | ✅ (subclase) |
| one-pol greedy 12 shells C-R-0014 | 9.313 Shor / 4.076 SOS | ✅ (subclase) |
| greedy 41 shells C-R-0012 | 9.504 Shor | ✅ (subclase) |
| **full dealias all-IC Sym** | **25.926** | ❌ |
| rank-1 Rayleigh (lower bound) | ≈3.45 | cota inf., no certifica |

**Conclusión:** barrera = gap de relajación Sym, no física. Stretch real ≪ C_† empíricamente.

---

## G. Escalera SOS / Shor (ataque C_†) — resultados

### G.1 Hermitian bands (2 polarizaciones, shells consecutivos)

| Shells | D | C_fullsym | C_ub SOS | Solver | vs C_† |
|--------|---|-----------|----------|--------|--------|
| {1,2,3} | 52 | 1.37 | **0.844** | Clarabel | ✅ |
| {1..4} | 64 | 1.37 | 0.894 | Clarabel | ✅ |
| {1..5} | 112 | 2.15 | 1.305 | Clarabel | ✅ |
| {1..6} | 160 | 2.99 | 1.701 | Clarabel | ✅ |
| {1..8} | 184 | 3.29 | 1.852 | COSMO | ✅ |
| {1..9} | 208 | — | 2.223 | COSMO | ✅ |
| {1..10} | 232 | — | 2.440 | COSMO | ✅ |
| {1..12} | 280 | — | 2.783 | COSMO | ✅ |
| **extrap full dealias** | — | 25.93 | **~14.7** | extrap | ❌ |

### G.2 One-pol (1 polarización/modo) — mejor subclase estructurada

| Witness | Shells | D | C_Shor | C_ub SOS | Status |
|---------|--------|---|--------|----------|--------|
| C-R-0008 | {1,2,3,4,5,6,8} | 92 | 8.476 | 3.356 | proved |
| C-R-0013 | greedy 11 shells | 154 | 9.189 | 4.015 | proved |
| C-R-0014 | +shell 24 (12 total) | 178 | 9.313 | 4.076 | proved |
| **techo +shell 25** | 13 shells | ~190 | **≈9.79** | — | **> C_†** |

**Techo one-pol:** shell 25 rompe C_†. No hay crecimiento greedy all-IC.

### G.3 Greedy monolítico (41 shells, C-R-0012)

- D ≈ 1772, C_Shor ≈ 9.504 ≤ C_† (Shor, **ya probado**)
- SOS extrap C_ub ≈ 5.39 (L-0055) — cierra subclase, **no all-IC**
- Solve estimado cluster: **~113 h** — deferred; **no cerraría C-0008 all-IC** aunque termine

### G.4 Caminos abandonados / techo

| Método | Resultado | Lemma |
|--------|-----------|-------|
| Column-triangle | C≈80 | L-0042 |
| A/B hybrid | C≈16.2 | L-0042 |
| Shell-block / Ky-Fan | C≳12–14 | L-0043 |
| min(defect, Young) | empeora | L-0026 |
| Inflar N-bound | no válido como majorant | L-0026 |
| Brute band SOS → full | extrap ~14.7 > C_† | L-0051 |

---

## H. 14 subclases probadas C-R-0002 … C-R-0014

Todas prueban Ω(0.02)≤41.284 **en su dominio** vía L-0026 (high slab) + L-0027 (low slab, C≤C_† en subclase).

| ID | Dominio (resumen) | C stretch | Lemma |
|----|-------------------|-----------|-------|
| C-R-0001 | mono-radial 1 shell | — | L-0020/21 |
| C-R-0002 | Ω₀≥Ω★ (high slab) | — | L-0024/26 |
| C-R-0003 | mono-radial α_r=0 | — | Stokes-Duhamel |
| C-R-0004 | shells {1,2,3} | 9.17 | L-0038 |
| C-R-0005 | 14 shells sparse | 9.52 | L-0039 |
| C-R-0006 | {1..5} | 9.01 | L-0040 |
| C-R-0007 | 8 shells con r=6 | 9.41 | L-0041 |
| C-R-0008 | one-pol {1..8} | 8.48 | L-0043 |
| C-R-0009 | rank-1 band | 2.99 | L-0044 |
| C-R-0010 | streaming | 9.12 | L-0045 |
| C-R-0011 | physical tensor | 9.32 | L-0044 |
| C-R-0012 | greedy 41 shells | 9.50 | L-0047 |
| C-R-0013 | one-pol greedy 11 | 9.19 | L-0066 |
| C-R-0014 | one-pol +shell 24 | 9.31 | L-0068 |

**Limitación:** subclases no se acumulan a all-IC. Son ataques estructurados hacia C-0008.

---

## I. L-0024 (pieza central C-0007)

En N=24, shell K²=147 sin self-triads bajo dealias ⇒ proyección u_K tiene N(u_K)=0.

Descomposición u = u_K + u_<, defecto δ = K²E − Ω ≥ gap·E_<.

Majorante:
```
‖N(u)‖ ≤ a(E)√δ + b·δ
```

Identidades enstrofía/energía ⇒ ODE comparación acoplada (E,Ω). Peor caso Ω(T)≈41.743 at Ω₀=0.5.

**Monotone in Ω₀:** Ω(T) no crece al aumentar Ω₀ ⇒ worst at Ω₀=E₀.

---

## J. Certificados N5 (59 JSON en repo)

Principales para revisión:

| Certificado | Contenido |
|-------------|-----------|
| CERT-L0070-all-ic-C0007-N24 | Cierre all-IC M=41.743 |
| CERT-L0024-spectral-defect-C0005-N24 | ODE defecto espectral |
| CERT-L0026-c0007-techo-N24 | Case-split Ω★ |
| CERT-L0027-low-slab-cubic-C0007-N24 | C_†, ODE cúbico |
| CERT-L0047-greedy-fullsym-C0007-N24 | C-R-0012 Shor |
| CERT-L0052 … L-0062 | SOS bands Hermitian |
| CERT-L0053 | one-pol C-R-0008 |
| CERT-L0065, CERT-L0067 | one-pol greedy SOS |

---

## K. Infra / herramientas usadas

- **Python:** lemmas N7, certificados, tests pytest
- **SOS:** Julia TSSOS order=2, sphere equality; Clarabel (D≤160), COSMO (D≤~180, ~5–6 min/run)
- **Lean N8:** esqueleto SosGram (native_decide C_ub < C_dagger para certs N5); ODE/PDE deferred
- **Cluster:** greedy D=1772 SOS ~113h — export ready, deferred

---

## L. Preguntas concretas para ChatGPT

1. **¿Cómo cerrar el gap 0.459** entre M=41.284 y L-0024=41.743 sin probar C≤C_† all-IC?
2. **¿Existe majorant rank-1 exacto** explotando estructura de triadas (p+q=k), antisimetría, div-free, sin pasar por Sym?
3. **¿SOS de grado bajo** en el polinomio `C_†²‖z‖⁶ − (2√2)²(z·L(zzᵀ))²` es tratable en D~6748 o requiere factorización por shells irreducible?
4. **¿Las 14 subclases C-R sugieren un patron** (one-pol, shells específicos) generalizable, o es combinatoria sin cierre?
5. **¿Relajar N** (p.ej. N=16) da un piloto cerrable que escale, o el techo Sym empeora igual?
6. **¿Buscar refutación** (IC adversarial con Ω(T)>41.284) es más rentable que certificar?
7. **¿El target M=41.284** (1.011× Stokes floor) es natural, o deberíamos declarar C-0007=41.743 como resultado final?

---

## M. Índice de archivos en el repo (para adjuntar)

### Documentos prioritarios (copiar/adjuntar a ChatGPT)

1. `reports/PAQUETE_CHATGPT_REVIEW_COMPLETO.md` — **este archivo**
2. `reports/C0007_STRUCTURED_CLOSURE_MAP.md` — mapa de cierre
3. `reports/ONEPOL_SCAN_L0064.md` — scanner one-pol + techo shell 25
4. `reports/RESUMEN_C0007_PARA_ANALISIS.md` — resumen previo (pre-L-0070, aún útil en §F–G)

### JSON clave

5. `conjectures/proved_restricted/C-0007.json`
6. `conjectures/active/C-0008.json`
7. `conjectures/proved_restricted/L-0070.json`
8. `conjectures/proved_restricted/L-0026.json`
9. `conjectures/proved_restricted/L-0027.json`
10. `conjectures/proved_restricted/L-0048.json` (techo fullsym)
11. `conjectures/proved_restricted/C-R-0014.json` (mejor one-pol)
12. `tools/sos_julia/data/scale_results.json` (toda la escalera SOS)

### Certificados

13. `certificates/CERT-L0070-all-ic-C0007-N24.json`
14. `certificates/CERT-L0067-sos-onepol-greedy-second-r24-N24.json`

### Código (si ChatGPT puede leer repo)

- `python/ns_exploration/conjectures/l0024_spectral_defect.py` — ODE L-0024
- `python/ns_exploration/conjectures/l0027_low_slab_cubic.py` — C_†, ODE cúbico
- `python/ns_exploration/conjectures/l0070_c0007_all_ic_l0024.py` — cierre Option A
- `python/ns_exploration/conjectures/l0064_onepol_scanner.py` — scanner subclases

---

## N. Prompt sugerido para ChatGPT

```
Eres un revisor matemático/numerico. Te adjunto un paquete sobre cotas de enstrofía
para Navier-Stokes GALERKIN TRUNCADO (N≤24, T³, dealias 2/3). NO es Clay ni continuo.

Estado: C-0007 all-IC CERRADO a Ω(0.02)≤41.743 (L-0024). C-0008 SHARP ABIERTO:
Ω(0.02)≤41.284. El cuello de botella histórico era probar C≤9.562 para el stretch
rank-1 en la franja baja Ω₀≤18.87; Sym full dealias da C≈25.9; SOS extrap ≈14.7;
subclases one-pol llegan a C≈9.31 (12 shells) pero shell 25 rompe C_†.

Por favor:
1) Resume el estado en 10 líneas.
2) Evalúa viabilidad de cerrar C-0008 con técnicas NO ya listadas como fallidas.
3) Propón 3 direcciones concretas nuevas (con ecuaciones/relajaciones específicas).
4) Di si conviene seguir computando o pivotar el enunciado.
5) Señala errores lógicos obvios en la cadena L-0024→L-0070→C-0007.

[Pegar secciones B–L de este documento, o adjuntar el .md completo]
```

---

*Sistema finito Galerkin. No continuo. No Clay. Repo: bastardus2 / NS-MRL Route B.*
