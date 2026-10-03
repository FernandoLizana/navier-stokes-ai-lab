# C-0008 Handoff Pack

**Fecha:** 2026-08-03  
**Estado:** `exploring` — NO promover a `proved`  
**Target sharp:** M ≈ 41.283999509509286

## Mejor bound conocido

| Lemma | Método | Ω(T)^hi | I_term^hi |
|-------|--------|---------|-----------|
| **L-0073** | Route A' equal_12 cluster+best | **54.23** | **37.92** |
| L-0075 | triad_min_defect (all-IC sketch) | 51.0 | — |

**Gap vs M:** +12.9 (terminal) / +9.7 (low-slab sketch)

## Contenido del pack

| Carpeta | Contenido |
|---------|-----------|
| `reports/` | Todos los reportes C-0008 + MASTER_HANDOFF |
| `certificates/` | CERT-L0072, CERT-L0073 |
| `conjectures/` | C-0008, L-0072..L-0075 |
| `lemmas/` | L-0071..L-0075 markdown |
| `experiments/terminal_weighted/` | Manifests, summaries, checkpoints, logs clave |
| `scripts/` | reproduce, resilient runners, verify |
| `python/experiments/` | Sprint scripts C-0008 |
| `python/tests/` | test_l0072..l0075, test_phase_d |
| `python/ns_exploration/terminal_weighted/` | Módulo core |
| `canvas/` | Canvas interactivo master handoff |
| `docs/` | Prompt Phase D |

## Reproducir (Windows)

```powershell
Set-Location <repo-root>
powershell -NoProfile -ExecutionPolicy Bypass -File .\scripts\reproduce_c0008_certificate.ps1
python -m ns_exploration.experiments.finalize_route_a_ladder
python -m ns_exploration.experiments.sprint_c0008_hybrid_low_slab
```

## Rutas refutadas (no reintentar igual)

- L-0072 per-shell Phase E integral
- L-0074 particiones más finas (equal_24/48, fine_low)
- Route D terminal peor que A' en optimum
- One-pol SOS beyond shell 25 (techo L-0069)

## Próximos pasos prioritarios

1. **A1** — Probar all-IC C ≤ C_† ≈ 9.56 (L-0027 stretch)
2. **A2** — Rigor high-slab L-0024/L-0026 + low-slab
3. **A3** — Nueva estructura terminal (non-partition)

Ver `reports/C0008_MASTER_HANDOFF.md` para detalle completo.
