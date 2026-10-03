# C-0008 Terminal-Weighted Scan (N2 exploratory)

> **C-0008 sigue abierta.** Esto no es certificado N5.

## Configuración

- Band radii: `[1, 2, 3, 4, 5, 6, 7, 8]` (no full dealias all-IC)
- Grid step: `0.001` (21 nodos; spec pedía 0.00025 — refinamiento pendiente)
- Evidencia: **N2 float**

## Resultados band {1..8}

| Métrica | Valor |
|---------|-------|
| **max_grid_bound** | **3.287687** at t=**0.02** |
| trapezoid_estimate | 0.064840 |
| left_riemann | 0.064840 |
| right_riemann | 0.064840 |
| **threshold integral** | **1.299313** |
| **margin vs threshold** | **+1.234473** |
| uniform threshold | 64.965638 |
| Closes integral (band)? | **True** (exploratorio) |
| Closes uniform (band)? | **True** (exploratorio) |

## Full dealias (L-0048 regression, t=T)

| Métrica | Valor |
|---------|-------|
| C_term^ub(T) full dealias | **25.925922** (L-0048 frozen / sparse {1..6} validated) |
| Integral if constant C≈25.9 | ≈ **0.52** ≪ 1.299 |
| max vs uniform threshold | 25.9 ≪ 64.97 |

**Honestidad:** el scan ejecutado es sobre **band {1..8}** (D≈184). All-IC requiere
full dealias (D=6748) en cada t con cobertura intervalar (Phase D). El valor terminal
25.93 ya conocido no supera el umbral uniforme, pero **falta certificar C(t) en todo [0,T]**.

## Recomendación

Proceder a **Phase D** (cobertura temporal certificada) sobre **full dealias**, no solo band.
