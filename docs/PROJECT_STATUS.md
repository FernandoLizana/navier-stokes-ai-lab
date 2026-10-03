# Project status (verified vs reported)

Last updated during the public repositioning review (local checkout).  
**Precedence:** live code and certificates under active paths > repair evaluation reports > historical handoff packs in `artifacts/`. Snapshots that say `RUNNING` describe a state at that date, not a live process.

## Status table

| Component / result | Source | Implemented? | Verification in this review | Scope / limits | Status |
|--------------------|--------|--------------|-----------------------------|----------------|--------|
| Public README positioning | `README.md`, `README.es.md` | Yes | Read + cross-checked with C-0008 / repair docs | Presentation only | **checked** |
| Pseudospectral solver | `python/ns_exploration/spectral/solver.py` | Yes | Demo + `test_sprint01` intended | Truncated \(\mathbb{T}^3\), float, ≤ N2 | **checked** (demo ran) |
| Taylor–Green / ABC / random IC | `python/ns_exploration/initial_conditions/` | Yes | Demo TG; tests exist | Div-free by construction + Leray | **observed in code** |
| Public demo CLI | `python -m ns_exploration.demo` | Yes | Ran N=8, t_end=0.02; wall ≈ 0.1 s | Budget-capped; no GPU/API | **checked** |
| C-0008 conjecture record | `conjectures/active/C-0008.json` | Registry JSON | Read fields | Galerkin N≤24; `exploring`; no Clay claim | **observed in code** |
| Repair Route A I_hi ≈ 170.59 | `reports/C0008_REPAIR_EVALUATION.md`, cert JSON | Reported | Numbers cited from report/cert metadata | Finite repair baseline; does not close C-0008 | **historically reported** |
| Cluster sketch I_hi ≈ 38.52 | `reports/L0073R_CLUSTER_AUDIT.md` | Reported + audit | Not re-run full cluster build | Formal N7 triangle gap remains | **historically reported** / sketch |
| Historical L-0073 ~37.92 | Same reports | Superseded invalid | Documented as invalid under repair | Do not treat as rigorous upper bound | **replaced / invalid** |
| Rational n=24 regen “RUNNING” | `reports/C0008_CERT_REPAIR_STATUS.md` | Campaign | Status text as of report date | Not assumed running now | **historically reported** |
| Band repair CI workflow | `.github/workflows/c0008-repair-band.yml` | Present | Not executed on GitHub here | Needs gmpy2; heavier than demo | **pending** (remote) |
| Lean NSGalerkin | `lean/NSGalerkin/` | Finite lemmas | README scope read; `lake build` not run here | No continuum PDE | **observed in code** |
| Lean BasicDefinitions stub | `lean/NavierStokes/BasicDefinitions.lean` | Stub | Path noted in prior review | Not a formalization of NS | **observed in code** |
| Julia / SOS | `julia/` | Present | Not executed | Optional stack | **pending** |
| Large artifact `l0049_full_M.npz` | `reports/` | Binary | Size noted previously (~hundreds of MB) | Regenerable / LFS candidate | **observed** |
| `compute_pack/` | `compute_pack/` | Helper for heavy jobs | Scripts present | Portable compute; not required for demo | **observed in code** |
| `artifacts/C0008_HANDOFF_PACK/` | `artifacts/` | Snapshot | Distinguish from live `python/` sources | Handoff freeze | **snapshot** |

## Precedence for contradictions

1. **Live package** under `python/ns_exploration/` and root install metadata.
2. **Active conjecture / certificate paths** (`conjectures/active/`, `certificates/`) plus adversarial verifiers—PASS means only the checks implemented.
3. **Repair / audit reports** (`reports/C0008_*`, `reports/L0073R_*`) for narrative of invalidation and gaps.
4. **Handoff packs / compute_pack data** as dated snapshots.
5. Older sprint READMEs that still say “Sprint 0–1 under construction” are **stale relative to later campaigns**; the public README supersedes them for newcomers.

Do not edit certificate payloads solely to change wording if hashes bind content; use docs like this file.
