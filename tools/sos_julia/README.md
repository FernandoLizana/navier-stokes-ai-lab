# tools/sos_julia — sparse SOS sidecar for C-0007 (Phase 5)

Attack all-IC (or large bands) via **TSSOS + Clarabel** on the physical
cubic `A(z,z,z)` with `‖z‖²=1`. Target: `C_ub = 2√2·ub_raw < C_† ≈ 9.562`
(equivalently `ub_raw < β ≈ 3.3807`).

## Status (2026-07-20)

- Julia **1.12.6** via `juliaup` (winget).
- Two envs:
  - `tools/sos_julia` — SumOfSquares experiments (odd-degree Ideal API still flaky)
  - `tools/sos_julia/tssos_env` — **TSSOS 1.5.3 + Clarabel** (working)
- Scale (order-2 sphere; Clarabel OK ≤D≈160; **COSMO** preferred for D≳180):

  | radii | D | C_fullsym | C_ub SOS | solver | beats Shor? |
  |-------|---|-----------|----------|--------|-------------|
  | 1–3 | 52 | ≈1.370 | ≈0.844 | Clarabel | yes |
  | 1–4 | 64 | ≈1.370 | ≈0.894 | Clarabel | yes |
  | 1–5 | 112 | ≈2.148 | ≈1.305 | Clarabel | yes |
  | 1–6 | 160 | ≈2.993 | ≈1.701 | Clarabel | yes |
  | 1–8 | 184 | ≈3.288 | ≈1.852 | COSMO | yes |
  | 1–9 | 244 | ≈3.991 | ≈2.223 | COSMO | yes |
  | 1–10 | 292 | ≈4.260 | ≈2.440 | COSMO | yes |
  | 1–12 | 356 | ≈4.912 | ≈2.783 | COSMO | yes |
  | one-pol C-R-0008 | 92 | ≈8.475 | ≈3.356 | Clarabel | yes (N5 L-0053) |

- See `data/scale_results.json`. Extrapolation: need SOS/Shor ratio ≲0.37 on full dealias to beat `C_†`; Hermitian band ratio ≈0.56 ⇒ **brute band SOS techo** (L-0051). One-pol structured SOS ratio ≈0.40 on D=92.
- **Pivot (L-0051):** **Done:** one-pol SOS D=92 (`C_ub≈3.36`, N5 Gram `L-0053`); band `{1,2,3}` N5 Gram `L-0052` (`C_ub_hi≈0.844`).
- `pin_ram` default **off** (caused ReadOnlyMemoryError on D=500).
- CLI Clarabel: `run_tssos_band.jl <dir> [order] [TS] [merge] [loose] [dualize]`
- CLI COSMO: `run_tssos_cosmo.jl <dir> [order] [TS] [merge]`
- Odd cubic ⇒ one SDP. Shell r=7 adds no modes.

## Install

```powershell
winget install Julialang.Juliaup --accept-package-agreements
# refresh PATH, then:
cd tools/sos_julia/tssos_env
julia --project=. bootstrap.jl
```

## Pipeline

```powershell
# 1) Export cubic from Python
cd repo_root
$env:PYTHONPATH='python'
python -m ns_exploration.experiments.export_sos_cubic --dense --radii 1 2 3 --out tools/sos_julia/data/band_123

# 2) Run TSSOS order-2
cd tools/sos_julia
julia --project=tssos_env run_tssos_band.jl data/band_123 2

# 3) N5 Gram export + certificate (band {1,2,3})
julia --project=tssos_env export_tssos_gram.jl data/band_123
python -m ns_exploration.experiments.sprint02_l0052_gram_n5

# 4) One-pol C-R-0008 N5 Gram (D=92, ~10 min SDP + export)
julia --project=tssos_env export_tssos_gram.jl data/onepol_cr0008
python -m ns_exploration.experiments.sprint02_l0053_gram_n5
```

## Honesty

- Float SDP without exact / rationalized Gram = **N2 numerical evidence**, not a theorem.
- Finite Galerkin only. No continuum / Clay claims.
- Full dealias (D=6748) still open; scale bands gradually.
