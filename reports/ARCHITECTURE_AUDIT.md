# Architecture Audit — Sprint 0/1

**Date:** 2026-07-16  
**Repo state at audit:** empty (greenfield). No pre-existing solver code to preserve.

## Decision

Create the full NS-MRL directory tree at the workspace root (`bastardus2`), matching the prescribed architecture. Implement only the Sprint 1 critical path; leave C++/CUDA/dashboard/ML as empty scaffolds.

## What exists after Sprint 1

| Area | Status |
|------|--------|
| Docs (Clay, foundations, spaces, criteria, integrity, roadmap, …) | Present |
| Literature matrix (≥40 primary sources) | Present |
| Python exploratory pseudospectral solver (T³) | Present + tested |
| Initial conditions (TG, ABC, random, tubes) | Present + tested |
| Python finite-mode interval probes | Present + tested |
| Julia ValidatedNS skeleton | Present (requires Julia) |
| Lean stubs | Present (no proofs) |
| C++/CUDA/dashboard | Empty dirs only |
| Pilot campaign | Runnable |

## Risks

- Dual Fourier conventions: mitigated by single module `fourier_conventions.py`.
- Treating float runs as proofs: blocked by evidence tags and integrity docs.
- Julia IntervalArithmetic not audited for FFT: Sprint 1 avoids FFT in Julia path.
