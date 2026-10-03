# C-0008 Compute Pack — portable CPU jobs (Ubuntu)

Todo el **estado** (checkpoints, logs, manifiestos, certs) vive en **`compute_pack/data/`**.  
Copia **`compute_pack/`** + el código **`python/`** (o el repo entero) al notebook; al terminar, trae solo **`compute_pack/data/`** de vuelta.

**C-0008 sigue `exploring`** — estos jobs endurecen certificados, no prueban cierre.

---

## 1. En tu PC (Windows) — exportar progreso actual

```powershell
cd C:\Users\...\bastardus2\compute_pack
.\export_state.ps1
```

Incluye checkpoints racional n=24 (**27/87** shells al último stop).

---

## 2. Qué copiar al notebook Ubuntu

| Opción | Qué copiar |
|--------|------------|
| **Mínimo** | `bastardus2/python/` + `bastardus2/compute_pack/` |
| **Recomendado** | Repo `bastardus2/` completo |

No hace falta GPU. Sí: **Python 3.11+**, **gmpy2** (MPFR).

### RAM alta (notebook con mucha RAM)

Workers se eligen solos según RAM libre + núcleos CPU:

```bash
source compute_pack/env.sh
./compute_pack/run.sh ram          # JSON: ram, cpu, recommended_workers
./compute_pack/run.sh rational_full_n24   # sin 2º arg → auto
```

| RAM libre | Workers típicos (~1.75 GB/worker) |
|-----------|-----------------------------------|
| 16 GB | ~8 |
| 32 GB | ~17 |
| 64 GB | ~32 (tope default `C0008_WORKERS_MAX`) |

Subir el tope si tienes 64 GB+ y muchos cores:

```bash
export C0008_WORKERS_MAX=48   # o 64
./compute_pack/run.sh rational_full_n24
```

Con **16 workers** el job P1 (~60 shells restantes) baja de **~16 h → ~6–8 h** reloj (si el CPU aguanta).

---

## 3. Setup Ubuntu (una vez)

```bash
cd bastardus2/compute_pack
chmod +x setup.sh run.sh env.sh pack_results.sh import_results.sh
./setup.sh
source env.sh
./run.sh list
```

---

## 4. Ejecutar jobs

### Un job (recomendado)

```bash
source compute_pack/env.sh
./compute_pack/run.sh rational_full_n24    # auto workers desde RAM
# o explícito:
./compute_pack/run.sh rational_full_n24 16
tail -f compute_pack/data/logs/supervisor_rational_full_n24.log
```

### Cola automática (prioridad 1→10, respeta dependencias)

```bash
./compute_pack/run.sh queue 8
```

### Estado / parar

```bash
./compute_pack/run.sh status
./compute_pack/run.sh stop rational_full_n24   # o stop sin id = todos
```

**Resiliencia:** cada job tiene checkpoint intra-shell, heartbeat, supervisor con restart si crashea o heartbeat >25 min, logs en `data/logs/`.

---

## 5. Traer resultados al PC principal

En el notebook:

```bash
./compute_pack/pack_results.sh
# → c0008_compute_results_YYYYMMDD.tar.gz
```

En Windows (PowerShell):

```powershell
cd bastardus2\compute_pack
bash import_results.sh c0008_compute_results_....tar.gz
```

---

## 6. Listado exhaustivo de jobs CPU

Tiempos = **horas reloj** con **auto workers** (ver `./run.sh ram`). Con **32 GB+ RAM** y CPU fuerte, P1 suele **~6–10 h** restantes (vs ~16 h con 3–4 workers).  
Ver también `./run.sh list`.

| Prioridad | Job ID | Qué hace | Horas | Deps | Resume |
|-----------|--------|----------|-------|------|--------|
| **1** | `rational_full_n24` | 87 shells full-dealias n=24, coeficientes MPFR racionales | **~6–16** (quedan ~16) | — | **27/87** |
| **2** | `cluster_rational_n24` | Manifiesto + cert L-0073R cluster desde manifest racional | **~6** | (1) | — |
| **3** | `rational_full_n20` | Full dealias n=20 racional (escala N) | **~10** | — | 0/87 |
| **4** | `rational_full_n16` | Full dealias n=16 racional | **~8** | — | 0/87 |
| **5** | `rational_full_n12` | Full dealias n=12 racional | **~6** | — | 0/87 |
| **6** | `float_full_n20` | Full dealias n=20 float (ladder L-0076, opcional) | **~10** | — | — |
| **7** | `float_full_n16` | Full dealias n=16 float | **~8** | — | — |
| **8** | `float_full_n12` | Full dealias n=12 float | **~6** | — | — |
| **9** | `cluster_float_n24_rebuild` | Rebuild cluster float (sanity) | **~5** | — | — |
| **10** | `band_rational_n24` | Band {1..6} racional (smoke, minutos) | **~0.1** | — | done |

### Total si corres TODO (secuencial)

| Paquete | Horas |
|---------|-------|
| Solo P1+P2 (n=24 racional + cluster) | **~28** |
| + escalas n=12,16,20 racional | **+24** → **~52** |
| + float n=12,16,20 opcional | **+24** → **~76** |

**Recomendación notebook:** P1 `rational_full_n24` → P2 `cluster_rational_n24`. El resto solo si sobra CPU.

---

## 7. Archivos en `compute_pack/data/`

```
data/
  terminal_weighted/          # checkpoints, manifests, shell_one_inf_ck/
  logs/                       # supervisor_*.log, child_*.log, nohup_*.log
  certificates/               # certs generados post-job
  export_manifest.json        # metadata export Windows
```

---

## 8. Variables de entorno (`env.sh`)

| Variable | Ruta |
|----------|------|
| `C0008_DATA_DIR` | `compute_pack/data/terminal_weighted` |
| `C0008_LOG_DIR` | `compute_pack/data/logs` |
| `C0008_CERT_DIR` | `compute_pack/data/certificates` |

El código en `python/` respeta `C0008_DATA_DIR` vía `portable_paths.py`.

---

## 9. Orden sugerido en el notebook

1. `./setup.sh`
2. `./run.sh rational_full_n24` — dejar overnight (auto RAM)
3. `./run.sh cluster_rational_n24 6`
4. `./pack_results.sh` — copiar tarball home

---

## 10. Qué NO esperar

- **No cierra C-0008** (I_hi target ~1.3; Route A ~170).
- Regen racional ≈ mismo I_hi que float (endurece N4, no magia numérica).
- Cluster racional requiere manifest 87/87 completo.
