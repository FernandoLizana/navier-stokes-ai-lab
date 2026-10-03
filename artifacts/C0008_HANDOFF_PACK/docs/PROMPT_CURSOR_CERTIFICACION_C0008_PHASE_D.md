# C-0008 — Fase D de certificación temporal full-dealias

## Especificación ejecutable para Cursor

Lee este documento completo antes de modificar código. Trátalo como una especificación de investigación matemática asistida por computador.

---

## 1. Estado de entrada

La primera entrega de la ruta terminal-weighted ya produjo:

```text
python/ns_exploration/terminal_weighted/
lemmas/L-0071-terminal-weighted-identity.md
reports/C0008_L0048_FILE_MAP.md
reports/C0008_TERMINAL_WEIGHTED_AUDIT.md
reports/C0008_TERMINAL_WEIGHTED_SCAN.md
experiments/terminal_weighted/repo_snapshot.json
experiments/terminal_weighted/scan.json
experiments/terminal_weighted/delivery_summary.json
python/ns_exploration/tests/test_terminal_weighted.py
python/ns_exploration/experiments/sprint_c0008_terminal_weighted_pilot.py
```

Resultados informados:

```text
Tests: 6 passed

Stokes floor hi:
40.824623079...

Full-dealias en t=T:
C_fullsym(T) ≈ 25.925922

Band {1..8}:
max_grid_bound ≈ 3.288 en t=T

Integral exploratoria band:
≈ 0.0648

Estimación heurística full-dealias suponiendo C constante:
≈ 0.52
```

La condición suficiente propuesta para cerrar C-0008 es:

\[
I_{\mathrm{term}}
:=
\int_0^T C_{\mathrm{term}}(t)\,dt
<
I_\star,
\]

con

```text
nu       = 0.1
E0       = 0.5
T        = 0.02
M_target = 41.283999509509286
```

y

\[
C_{\mathrm{term}}(t)
=
2\sqrt 2
\sup_{\|z\|_2=1}
\left|
\left\langle
A e^{-2\nu A(T-t)}z,\ N(z)
\right\rangle
\right|.
\]

El umbral nominal informado es:

```text
I_star ≈ 1.29931275566075
```

No copiar este decimal como axioma. Recalcularlo con aritmética rigurosa desde:

\[
I_\star
=
2\sqrt2
\left(
M_{\mathrm{target}}-\mathrm{StokesFloor}_{\mathrm{hi}}
\right).
\]

C-0008 continúa abierta hasta completar una cobertura certificada full-dealias y un verificador independiente.

---

# 2. Objetivo de esta fase

Construir una prueba reproducible de una de estas dos afirmaciones:

## Éxito principal

\[
I_{\mathrm{term}}^{\mathrm{hi}}
<
I_\star^{\mathrm{lo}},
\]

y por consiguiente:

\[
\Omega(0.02)
\le
\mathrm{StokesFloor}_{\mathrm{hi}}
+
\frac{I_{\mathrm{term}}^{\mathrm{hi}}}{2\sqrt2}
<
41.283999509509286.
\]

## Resultado honesto si no cierra

Entregar la mejor cota rigurosa obtenida, identificar exactamente qué componente impide cerrar y no cambiar el estado de C-0008.

---

# 3. Prohibiciones

No hacer ninguna de estas cosas:

1. No usar el scan de la banda `{1..8}` como prueba all-IC.
2. No asumir que `Cterm(t)` alcanza su máximo en `t=T`.
3. No asumir que `Cterm(t) <= Cterm(T)`.
4. No usar `25.925922 * T` como cota rigurosa sin demostrar uniformidad.
5. No usar trapecio float como certificado.
6. No describir una norma calculada por iteración de potencia como cota superior si solo entrega una cota inferior.
7. No llamar “certificado” a una salida SDP o SVD float sin control de error.
8. No usar una matriz full-dealias distinta de la que corresponde al lema L-0071.
9. No afirmar que probar `N=24` cubre automáticamente todo `N<=24` sin demostrar la relación entre truncaciones.
10. No modificar `C-0008.json` a `proved` antes de que el verificador independiente devuelva `PASS`.
11. No afirmar implicaciones para el PDE continuo ni para el problema Clay.
12. No ocultar redondeos, tolerancias, fallos de convergencia o dependencias opcionales.

---

# 4. Auditoría matemática previa obligatoria

Antes de calcular la integral, verificar por escrito y mediante tests:

## 4.1 Identidad terminal

Partiendo de:

\[
\dot u=-\nu Au-N(u),
\]

\[
W(t)=Ae^{-2\nu A(T-t)},
\]

\[
\Phi(t)=\frac12\langle u,W(t)u\rangle,
\]

demostrar:

\[
\Phi'(t)
=
-\langle W(t)u,N(u)\rangle,
\]

y:

\[
\Phi(T)=\Omega(T).
\]

Comprobar cuidadosamente:

```text
signo de N
factor 1/2
producto real/complejo
normalización Fourier
volumen de T^3
definición de A
definición de Omega
realidad u_-k = conjugado(u_k)
proyección de Leray
```

## 4.2 Stokes floor

Recalcular rigurosamente:

\[
\Phi(0)
\le
E_0
\max_{r\in\mathcal R_N}
r e^{-2\nu rT}.
\]

Entregar:

```text
shell r que maximiza
intervalo certificado del máximo
lista exacta de shells consideradas
prueba de que no falta ninguna shell
```

## 4.3 Relación entre fullsym y la forma cúbica

Verificar que el objeto acotado por la rutina fullsym satisface exactamente:

\[
| \langle W(t)z,N(z)\rangle |
\le
\frac{C_{\mathrm{term}}^{\mathrm{hi}}(t)}{2\sqrt2}
\|z\|_2^3.
\]

No aceptar una constante que corresponda a otra base, otra banda, otra polarización o una normalización distinta.

## 4.4 Cobertura de `N <= 24`

Determinar cuál de estas opciones es correcta:

1. una prueba para el espacio `N=24` domina todas las truncaciones menores;
2. las truncaciones no son anidadas bajo la convención usada;
3. hay que certificar cada `N` admisible por separado.

Si se usa la opción 1, demostrarla formalmente, incluyendo el efecto de la proyección Galerkin y el dealias 2/3.

---

# 5. Representación offline/online por shells

La dependencia temporal solo aparece mediante los pesos:

\[
w_r(t)=r e^{-2\nu r(T-t)}.
\]

Construir una descomposición exacta:

\[
\mathcal M(t)
=
\sum_{r\in\mathcal R}
\varphi_r(t)\mathcal M_r,
\]

donde:

```text
phi_r(t) = exp(-2*nu*r*(T-t))
```

y `M_r` contiene la contribución completa con el factor `r`, o usar una convención equivalente documentada.

Requisitos:

1. `M(t)` debe ser el mismo operador fullsym usado para `Cterm`.
2. Reconstruir `M(T)` desde los bloques y recuperar L-0048.
3. Verificar la reconstrucción en varias bandas pequeñas de forma densa.
4. Guardar hashes de cada bloque.
5. No almacenar una matriz densa imposible si existe una representación streaming.
6. Registrar dimensión de dominio y codominio.
7. Registrar número de shells y de coeficientes no nulos.

Crear:

```text
python/ns_exploration/terminal_weighted/shell_decomposition.py
experiments/terminal_weighted/shell_manifest.json
reports/C0008_SHELL_DECOMPOSITION_AUDIT.md
```

---

# 6. Derivadas temporales

Usar:

\[
w_r'(t)
=
2\nu r^2e^{-2\nu r(T-t)},
\]

\[
w_r''(t)
=
4\nu^2r^3e^{-2\nu r(T-t)}.
\]

Construir operadores:

\[
\mathcal M'(t),
\qquad
\mathcal M''(t).
\]

Agregar tests por diferencias finitas únicamente como regresión N2, nunca como prueba.

La certificación debe usar las fórmulas analíticas de los pesos.

---

# 7. Estrategias de certificación

Implementar las rutas en el orden siguiente. Cada ruta debe poder activarse mediante configuración y producir un manifiesto independiente.

## Ruta A — Cota integrada por shells

Si existen cotas certificadas:

\[
\|\mathcal M_r\|_{\mathrm{op}}\le B_r,
\]

entonces:

\[
\|\mathcal M(t)\|_{\mathrm{op}}
\le
\sum_r e^{-2\nu r(T-t)}B_r.
\]

Integrar exactamente:

\[
\int_0^T C_{\mathrm{term}}(t)\,dt
\le
2\sqrt2
\sum_r
B_r
\frac{1-e^{-2\nu rT}}{2\nu r}.
\]

Ajustar la fórmula si `r` ya está incluido o no en `M_r`. Documentar la convención y añadir una prueba dimensional.

Ventaja: evita discretización temporal.

Riesgo: la desigualdad triangular entre shells puede ser demasiado laxa.

Crear:

```text
python/ns_exploration/terminal_weighted/certify_shell_integral.py
experiments/terminal_weighted/shell_integral_certificate.json
reports/C0008_SHELL_INTEGRAL_CERTIFICATE.md
```

## Ruta B — Centros temporales + Lipschitz certificado

Para cada celda:

\[
I_j=[a_j,b_j],
\qquad
m_j=\frac{a_j+b_j}{2},
\qquad
h_j=b_j-a_j,
\]

obtener:

\[
\|\mathcal M(m_j)\|_{\mathrm{op}}
\le U_j.
\]

Obtener además una cota:

\[
\sup_{t\in I_j}
\|\mathcal M'(t)\|_{\mathrm{op}}
\le L_j.
\]

Entonces:

\[
C_{\mathrm{term}}(t)
\le
2\sqrt2
\left(
U_j+|t-m_j|L_j
\right),
\]

y:

\[
\int_{I_j}C_{\mathrm{term}}(t)\,dt
\le
2\sqrt2
\left(
h_jU_j+\frac{h_j^2}{4}L_j
\right).
\]

La forma más segura de obtener `L_j` es:

\[
L_j
\le
\sum_r
\sup_{t\in I_j}|w_r'(t)|
B_r,
\]

con `B_r` certificados.

No asumir monotonía de la norma total.

## Ruta C — Taylor de segundo orden

Si la Ruta B es demasiado laxa, usar:

\[
\mathcal M(t)
=
\mathcal M(m_j)
+
(t-m_j)\mathcal M'(m_j)
+
R_2(t),
\]

con:

\[
\|R_2(t)\|
\le
\frac{|t-m_j|^2}{2}
\sup_{s\in I_j}\|\mathcal M''(s)\|.
\]

Esto permite:

\[
\|\mathcal M(t)\|
\le
U_j
+
|t-m_j|U_j'
+
\frac{|t-m_j|^2}{2}Q_j.
\]

Integrar analíticamente cada término.

## Ruta D — Perturbación desde `t=T`

Usar:

\[
\mathcal M(t)
=
\mathcal M(T)+\Delta\mathcal M(t),
\]

\[
\|\mathcal M(t)\|
\le
\|\mathcal M(T)\|
+
\|\Delta\mathcal M(t)\|.
\]

Acotar el segundo término mediante los bloques por shell:

\[
\|\Delta\mathcal M(t)\|
\le
\sum_r
\left|
e^{-2\nu r(T-t)}-1
\right|B_r.
\]

Integrar esta desigualdad exactamente.

Esta ruta puede ser útil porque existe una cota terminal favorable, pero no debe confundirse con una demostración de monotonía.

---

# 8. Cota certificada de la norma de operador

Esta es la parte más delicada.

## 8.1 Auditoría de L-0048

Determinar si `25.925922...` es:

```text
a) valor aproximado de la norma;
b) cota inferior;
c) cota superior rigurosa;
d) intervalo certificado;
e) una cota analítica más laxa.
```

No continuar llamándolo “fullsym certificado” hasta documentar el mecanismo.

## 8.2 Requisito

Cada `U_j`, `B_r`, `L_j` o `Q_j` usado en la prueba final debe ser una **cota superior** con redondeo dirigido.

Métodos admisibles:

- aritmética de intervalos/bolas;
- factorización intervalar;
- cota residual rigurosa con hipótesis verificadas;
- norma `sqrt(||M||_1 ||M||_infinity)` calculada hacia arriba;
- Frobenius calculada hacia arriba;
- Lanczos verificado con una cota superior matemáticamente justificada;
- certificado racional;
- combinación de las anteriores.

Métodos no admisibles por sí solos:

- `numpy.linalg.svd`;
- `scipy.sparse.linalg.svds`;
- iteración de potencia;
- randomized SVD;
- un autovalor float sin residual riguroso;
- tolerancia del solver tomada como prueba.

## 8.3 Backend aritmético

Detectar qué backend riguroso está disponible.

Preferencias:

1. Arb mediante `python-flint`;
2. MPFR con redondeo dirigido mediante `gmpy2`;
3. Julia `IntervalArithmetic.jl`;
4. racionales exactos cuando sea viable.

No usar `decimal` o `mpmath` como si garantizaran automáticamente intervalos dirigidos.

Registrar:

```text
backend
version
precision
rounding mode
platform
```

---

# 9. Malla adaptativa

Crear un certificador adaptativo.

Algoritmo:

```text
queue = [[0,T]]
total_hi = 0

while queue not empty:
    I = pop(queue)
    compute rigorous integral upper bound on I

    if local bound quality is acceptable:
        accept I
    else:
        bisect I

stop when:
    total_hi < I_star_lo
or:
    resource limit reached
or:
    minimum interval width reached
```

No hace falta minimizar la integral con gran precisión. Basta cerrar con margen verificable.

Priorizar celdas donde:

```text
local integral contribution is large
derivative bound is loose
operator-norm uncertainty is large
```

Configuración sugerida:

```text
initial_cells = 4
max_cells = configurable
min_width = configurable
precision_bits = configurable
safety_margin = configurable
```

Crear:

```text
python/ns_exploration/terminal_weighted/adaptive_time_certificate.py
python/ns_exploration/experiments/sprint_c0008_phase_d.py
experiments/terminal_weighted/phase_d_config.json
```

---

# 10. Verificador independiente

Crear:

```text
verify_c0008_terminal_certificate.py
```

El verificador no debe confiar en el script que generó el certificado.

Debe:

1. cargar el manifiesto;
2. verificar hashes;
3. verificar parámetros matemáticos;
4. verificar cobertura exacta de `[0,T]`;
5. detectar huecos o solapamientos inválidos;
6. recalcular `I_star`;
7. recalcular Stokes floor;
8. verificar cada cota local;
9. sumar hacia arriba las contribuciones;
10. verificar:

\[
I_{\mathrm{term}}^{\mathrm{hi}}
<
I_\star^{\mathrm{lo}};
\]

11. verificar la cota final de enstrofía;
12. imprimir `PASS` o `FAIL`;
13. devolver código de salida distinto de cero al fallar.

No aceptar un manifiesto que incluya directamente los resultados finales sin datos suficientes para recomputarlos.

---

# 11. Formato del certificado

Crear:

```text
certificates/CERT-L0072-C0008-terminal-full-dealias.json
```

Campos mínimos:

```json
{
  "claim": "C-0008",
  "status": "candidate_certificate",
  "parameters": {
    "nu": "...",
    "E0": "...",
    "T": "...",
    "N_max": 24,
    "M_target": "..."
  },
  "normalization": {},
  "mode_set": {},
  "arithmetic_backend": {},
  "stokes_floor_interval": ["lo", "hi"],
  "i_star_interval": ["lo", "hi"],
  "method": "shell_integral | adaptive_lipschitz | second_order | terminal_perturbation",
  "operator_bound_method": "...",
  "time_cells": [],
  "integral_interval": ["lo", "hi"],
  "omega_T_interval": ["lo", "hi"],
  "strict_margin": "...",
  "source_hashes": {},
  "data_hashes": {},
  "verification_command": "...",
  "limitations": []
}
```

Solo después de `PASS`, generar una copia con:

```json
"status": "verified"
```

---

# 12. Tests obligatorios

Agregar al menos:

## Algebra

- identidad terminal simbólica;
- derivadas de pesos;
- reconstrucción por shells;
- igualdad band dense vs streaming;
- recuperación de L-0048 en `t=T`.

## Cobertura

- unión de celdas igual a `[0,T]`;
- sin huecos;
- orden correcto;
- suma intervalar hacia arriba.

## Seguridad

- modificar un hash produce `FAIL`;
- reducir artificialmente una cota produce `FAIL`;
- eliminar una celda produce `FAIL`;
- cambiar `nu`, `T`, `E0` o `M_target` produce `FAIL`;
- usar un resultado band como full-dealias produce `FAIL`.

## Regresión pequeña

Para una banda de dimensión manejable:

1. construir matriz densa exacta o de alta precisión;
2. comparar la cota certificada con SVD de alta precisión;
3. comprobar que la cota superior realmente queda por encima;
4. probar varios tiempos.

---

# 13. Criterios de éxito

## Certificación completa

Declarar C-0008 demostrada para el sistema truncado únicamente si:

```text
all mandatory tests pass
full-dealias is covered
all required N<=24 are covered
interval arithmetic is rigorous
independent verifier prints PASS
integral_hi < i_star_lo
omega_hi < M_target
strict margin > 0
```

## Éxito parcial

Si no cierra, entregar:

```text
best integral upper bound
best final Omega upper bound
largest local contribution
largest source of overestimation
required improvement factor
runtime and memory
recommended next method
```

## Refutación de la ruta, no de C-0008

Si una cota inferior rigurosa del mejor majorante elegido supera el umbral, solo concluir:

```text
this certification route cannot close C-0008 at the current relaxation level
```

No concluir que C-0008 es falsa.

---

# 14. Entregables finales

Crear:

```text
lemmas/L-0072-terminal-full-dealias-bound.md
reports/C0008_PHASE_D_CERTIFICATION.md
reports/C0008_PHASE_D_FAILURE_ANALYSIS.md        # solo si falla
certificates/CERT-L0072-C0008-terminal-full-dealias.json
verify_c0008_terminal_certificate.py
reproduce_c0008_certificate.ps1
reproduce_c0008_certificate.sh
experiments/terminal_weighted/phase_d_run.json
```

El informe debe responder:

```text
A. ¿Qué cantidad se certificó exactamente?
B. ¿La cota es full-dealias?
C. ¿Qué valores de N están cubiertos?
D. ¿Qué backend riguroso se usó?
E. ¿Cómo se acotó la norma de operador?
F. ¿Cuál es Iterm_hi?
G. ¿Cuál es Istar_lo?
H. ¿Cuál es el margen?
I. ¿Cuál es Omega(T)_hi?
J. ¿El verificador independiente devuelve PASS?
K. ¿Qué parte sigue siendo N2, N5, N7?
L. ¿Cuáles son las limitaciones?
```

---

# 15. Orden de ejecución

Trabaja en este orden:

```text
1. auditar L-0048
2. confirmar normalización
3. confirmar cobertura N<=24
4. construir bloques por shell
5. reproducir t=T
6. obtener B_r certificados
7. intentar Ruta A
8. si no cierra, Ruta D
9. si no cierra, Ruta B adaptativa
10. si no cierra, Ruta C
11. generar certificado
12. ejecutar verificador independiente
13. solo entonces actualizar el estado de C-0008
```

Debido al margen heurístico favorable, no ejecutar primero un SOS de mayor orden.

---

# 16. Mensaje de inicio que debes devolver

Antes de implementar, responde con:

1. clasificación exacta del valor L-0048;
2. método propuesto para cotas superiores rigurosas;
3. prueba o rechazo de que `N=24` domina `N<24`;
4. estimación de memoria para la descomposición por shells;
5. ruta A/B/C/D que intentarás primero;
6. archivos que modificarás;
7. riesgos matemáticos detectados.

Después comienza la implementación sin declarar C-0008 cerrada.
