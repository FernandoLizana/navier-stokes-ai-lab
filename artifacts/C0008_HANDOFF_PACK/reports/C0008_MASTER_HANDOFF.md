# C-0008 — Artefacto maestro (handoff completo)

**Fecha:** 2026-08-03  
**Estado conjectura:** `exploring` — **NO promover a `proved`**  
**Dominio:** Galerkin T³, N≤24, dealias 2/3 — finito-dimensional, no Clay

---

## 1. Objetivo sharp

| Quantity | Target | Stokes floor |
|----------|--------|--------------|
| **M** = Ω(T) bound | **41.283999509509286** | 40.82462307930434 |
| **I_*** (integral terminal) | < **1.2993** | — |

C-0007 ya probado en M≈41.743. C-0008 pide el refinamiento sharp (mid floor / L-0024).

---

## 2. Mapa de rutas (resultado)

```
                    C-0008 sharp M ≈ 41.28
                           │
         ┌─────────────────┼─────────────────┐
         │                 │                 │
    TERMINAL           LOW-SLAB          SOS / TRIAD
    (weighted)         (L-0027)          (subclass)
         │                 │                 │
    L-0072 REFUTED    L-0075 sketch      14 C-R proved
    per-shell         Ω≈51 all-IC        C_ub≈0.84 band
         │                 │                 │
    L-0073 BEST       gap +9.7           extrap ~17.5
    Ω≈54.2            vs M               (not all-IC)
         │
    L-0074 REFUTED
    finer clusters
```

---

## 3. Cronología de compute (realizado)

| Fase | Qué | Tiempo | Resultado |
|------|-----|--------|-----------|
| **Phase D** | Frobenius full-dealias 87 shells | ~39 min | I_hi≈112.8 Route D |
| **Phase E** | Per-shell min(F,1-inf) 87/87 | ~24h+ | I_hi≈78.8, Ω≈68.7 |
| **Route A'** | 12 clusters + best | ~69 min | **I_hi≈37.9, Ω≈54.2** (L-0073) |
| **Route A ladder** | equal_12/24/48 + fine_low | ~5h total | equal_12 gana (L-0074) |
| **Hybrid low-slab** | triad/defect/SOS ODE | ~1 min | all-IC Ω≈51.0 (L-0075) |
| **SOS/triad audit** | cert scan + weighted | ~2 min | all-IC open |

---

## 4. Mejores bounds por ruta

### 4.1 Terminal (full dealias, N5)

| Lemma | Método | I_term^hi | Ω(T)^hi | Verifier |
|-------|--------|-----------|---------|----------|
| L-0072 | Phase E Route D per-shell | 78.80 | 68.69 | FAIL |
| **L-0073** | **Route A' equal_12 cluster+best** | **37.92** | **54.23** | **FAIL** |
| L-0074 | equal_24 | 63.85 | 63.40 | — |
| L-0074 | fine_low | 40.13 | 55.01 | — |
| L-0074 | equal_48 (Route D) | 70.50 | 65.75 | — |

**Gap vs M:** Ω ≈ 54.2 − 41.28 ≈ **+12.9** (~31% above target)  
**Gap vs I_*:** I ≈ 37.9 vs 1.30 ≈ **29×**

### 4.2 Route A ladder full (2026-08-03, recomputado)

| Config | clusters | Route A | Route D | pick | Ω | runtime |
|--------|----------|---------|---------|------|---|---------|
| **equal_12** | 11 | 37.92 | 62.48 | A | **54.23** | ~72 min |
| fine_low | 19 | 40.13 | 62.27 | A | 55.01 | ~85 min |
| equal_24 | 22 | 63.85 | 65.58 | A | 63.40 | ~75 min |
| equal_48 | 44 | 103.87 | 70.50 | D | 65.75 | ~72 min |

**Conclusión L-0074:** 12 clusters equal es óptimo; más clusters empeoran.

### 4.3 Low-slab hybrid (L-0075)

| Ruta | Ω(T) worst (Ω₀≤Ω★) | All-IC? |
|------|---------------------|---------|
| triad_hybrid ODE | 144.16 | sketch |
| **triad_min_defect ODE** | **51.00** | **sketch all-IC** |
| weighted_min_defect | 57.30 | sketch |
| sos_greedy one-pol cubic | 21.71 | **subclase** |
| sos_band123 cubic | 17.54 | **subclase** |

- **All-IC sketch gap:** +9.7 vs M (mejor que terminal 54.2 en low slab solo)
- **Subclass SOS:** Ω≈17.5 — **no cierra C-0008 all-IC**

### 4.4 Stretch / SOS (all-IC open)

| Quantity | Value | vs C_†≈9.56 |
|----------|-------|-------------|
| C_emp (N2) | 0.0093 | yes (empirical) |
| C_fullsym all-IC | 25.93 | **no** (techo) |
| Best SOS subclass | 0.844 | subclass |
| Greedy one-pol SOS | 4.076 | subclass (C-R-0014) |
| SOS extrap D=6748 | ~17.5 | extrap, not proved |
| Weighted triad C | ~28.4 | no |

---

## 5. Lemmas y certificados

| ID | Status | Rol |
|----|--------|-----|
| L-0072 | route_refuted | Per-shell triangle + Phase E no cierra |
| **L-0073** | **best_known_bound** | Terminal cluster+best equal_12 |
| L-0074 | partition_ladder_refuted | Finer clusters no ayudan |
| L-0075 | best_low_slab_hybrid | Low-slab min(triad,defect,SOS) sketch |

### Certificados clave

- `certificates/CERT-L0072-C0008-terminal-full-dealias.json`
- `certificates/CERT-L0073-C0008-route-a-prime-best.json`

### Manifests clave

- `experiments/terminal_weighted/shell_manifest_best_full.json` (Phase E, 87 shells)
- `experiments/terminal_weighted/shell_manifest_cluster_best_full.json` (L-0073)
- `experiments/terminal_weighted/route_a_manifest_equal_{12,24,48}.json`
- `experiments/terminal_weighted/route_a_manifest_fine_low.json`
- `experiments/terminal_weighted/hybrid_low_slab_summary.json`

---

## 6. Rutas REFUTADAS (no reintentar igual)

1. **L-0072** — per-shell triangle + min(F,1-inf) integral full-dealias
2. **L-0074** — particiones adaptativas más finas (24/48/fine_low)
3. **Terminal Route D cluster** — peor que A' en optimum
4. **One-pol SOS beyond shell 25** — techo L-0069 (Shor > C_†)

---

## 7. Scripts resilientes y checkpoint

| Script | Checkpoint | Resume |
|--------|------------|--------|
| `run_upgrade_best_resilient.ps1` | `upgrade_best_checkpoint.json` + per-shell ck | auto |
| `run_route_a_ladder_resilient.cmd` | `route_a_ladder_full_checkpoint.json` | `--resume` |
| `reproduce_c0008_certificate.ps1` | — | verify-only ~3 min |

**Nota:** pase `cluster best` D² **no** tiene checkpoint intra-chunk (~70 min por config entera).

### Comandos reproduce

```powershell
Set-Location "<REPO_ROOT>"
powershell -NoProfile -ExecutionPolicy Bypass -File .\reproduce_c0008_certificate.ps1
python -m ns_exploration.experiments.finalize_route_a_ladder
python -m ns_exploration.experiments.sprint_c0008_hybrid_low_slab
```

---

## 8. Posibles soluciones (priorizadas)

### A. Alta prioridad / estructura nueva

| # | Ruta | Idea | Esfuerzo | Cierra all-IC? |
|---|------|------|----------|----------------|
| **A1** | **L-0027 all-IC** | Probar C ≤ C_†≈9.56 en low slab + L-0026 high slab | Matemático + empírico | **Sí, si probado** |
| **A2** | **High-slab L-0024/L-0026** | Combinar slabs rigurosamente para Ω≤M | Medio | Parcial |
| **A3** | **Nueva estructura terminal** | No-partition blocks, new inequality class | Alto | Posible |

### B. Media prioridad / compute

| # | Ruta | Idea | Tiempo est. | Nota |
|---|------|------|-------------|------|
| B1 | Weighted triad refinement | Mejor peso que inv_sqrt_rp_rq | Horas | C~28 hoy |
| B2 | SOS band extend | band12345678910+ extrap riguroso | Días | subclass |
| B3 | N-ladder full | n=12,16,20,24 full-dealias domination | Días×4 | pilot band hecho |
| B4 | Checkpoint intra-chunk cluster | Evitar recomputar 70 min | Dev 1 día | infra |

### C. Baja prioridad / refutado o marginal

| # | Ruta | Por qué no |
|---|------|------------|
| C1 | Más clusters equal_24/48 | L-0074 refutado |
| C2 | Route D terminal | Peor que A' |
| C3 | One-pol SOS shell≥25 | Techo L-0069 |
| C4 | Cluster greedy SOS D=1772 | ~113h, subclass, extrap only |

### D. Infra / operación (PC viejito)

- Usar `*_resilient*` scripts, `--workers 1-2`
- Repo fuera de OneDrive si es posible
- Desactivar sleep: `powercfg /change standby-timeout-ac 0`
- PowerShell externo, no Cursor background
- Desactivar NVIDIA overlay

---

## 9. Gap analysis honesto

| Obstacle | Detail |
|----------|--------|
| Terminal integral | Best Ω≈54.2 >> M≈41.28 |
| I_* strict | Best I≈37.9 >> I_*≈1.30 |
| All-IC C | C_fullsym≈25.9 >> C_†≈9.56 |
| Low-slab only | Best all-IC sketch Ω≈51 still +9.7 |
| Subclass SOS | Looks good (Ω≈17.5) but **not all initial conditions** |

**C-0008 permanece abierto.** El pipeline terminal está agotado en la familia cluster/partition probada.

---

## 10. Tests

```powershell
Push-Location python
python -m pytest ns_exploration/tests/test_phase_d.py `
  ns_exploration/tests/test_l0073_route_a_prime.py `
  ns_exploration/tests/test_l0074_route_a_ladder.py `
  ns_exploration/tests/test_l0075_hybrid_low_slab.py -q
Pop-Location
```

**Esperado:** 21+ tests PASS; verifiers L-0072/L-0073 **FAIL** (correcto).

---

## 11. Índice de reportes

| Reporte | Contenido |
|---------|-----------|
| `C0008_PHASE_D_CERTIFICATION.md` | Phase D Frobenius |
| `C0008_ROUTE_A_PRIME_BEST.md` | L-0073 |
| `C0008_ROUTE_A_LADDER.md` / `_FULL.md` | Ladder |
| `C0008_HYBRID_LOW_SLAB.md` | L-0075 |
| `C0008_SOS_TRIAD_STRETCH.md` | SOS audit |
| `C0008_WEIGHTED_TRIAD_SOS.md` | L-0031 |
| `SPRINT_C0008_L0027_POST_TERMINAL.md` | Stretch audit |
| **`C0008_MASTER_HANDOFF.md`** | **Este documento** |

---

*Generado como handoff completo del sprint C-0008 terminal + ladder + hybrid. Verifier FAIL es el comportamiento esperado hasta cerrar M sharp.*
