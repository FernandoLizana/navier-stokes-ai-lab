# Roadmap

Evidence discipline: simulations ≤ N3; proofs for truncated Galerkin ≤ N7–N8;
continuum Clay claims only after a *general* technique exists. No Clay claim yet.

| Phase | Focus | Status |
|-------|-------|--------|
| 0 | Foundations, Clay target, literature | **Done** |
| 1 | Exploratory pseudospectral solver + tests | **Done** (Sprint 1) |
| 2 | Adversarial search, residuals, adjoint, CMA | **Done** (Sprint 2) |
| 3 | Finite Galerkin lemmas L-0001…L-0006 + N5 certificates | **Done** (partial Sprint 2–3) |
| 4 | Conjecture line C-000x / C-S-000x | **C-0004…7** proved; **C-0008** sharp open; **C-R-0002…0014** proved |
| 5 | Bridge lemmas (L-0007…L-0070) | **L-0070** all-IC C-0007 at M≈41.743; C-0008 sharp open |
| 6 | Lean N8 (structural + certificate arithmetic) | **Skeleton done**; ODE/PDE deferred |
| 7 | CAP / Newton–Kantorovich a posteriori | Deferred |
| 8 | Julia validated-numerics stack | Stub (Python N5 covers finite identities) |
| 9 | Partial publication of truncated results | Deferred |
| 10 | Clay-route evaluation (B primary, D secondary) | Only after general continuum technique |

## Active conjectures (N6, unproved)

- **C-0008** — sharp refinement of C-0007 (M≈41.284); **L-0056** closure map (`reports/C0007_STRUCTURED_CLOSURE_MAP.md`)

## Proved all-IC (N7)

- **C-0007** — all IC, Ω(0.02) ≤ **41.74337594** via **L-0070** / L-0024 spectral-defect (Option A)

## Proved restricted (N7) + certified constants (N5)

- L-0001…L-0006 (hand, truncated)
- L-0007…L-0070 (**L-0070** all-IC C-0007; **L-0069** laptop one-pol; **C-R-0014**; cluster deferred)
- CERT-L0006-shell-k2; CERT-L0003-shell-k3; CERT-L0003-full-dealias-N12/N16
- CERT-L0018-envelope-CS0002-N16; CERT-L0018-envelope-C0003-N16; CERT-L0018-envelope-C0004-N24
- CERT-L0019-stokes-floor-N24; CERT-L0020-duhamel-H1-N24; CERT-L0021-quartic-radial-N24
- CERT-L0022-gamma-bootstrap-N24; CERT-L0023-cubic-dissipation-N24
- CERT-L0024-spectral-defect-C0005-N24; CERT-L0025-spectral-defect-C0006-N32
- CERT-L0026-c0007-techo-N24; CERT-L0027-low-slab-cubic-C0007-N24
- CERT-L0070-all-ic-C0007-N24
- CERT-L0028-embedding-techo-C0007-N24; CERT-L0029-frobenius-stretch-C0007-N24
- CERT-L0030-triad-N-C0007-N24; CERT-L0031-0033-techos-C0007-N24
- CERT-L0034-shell-N-C0007-N24
- CERT-L0035-ginf-stretch-C0007-N24
- CERT-L0036-signed-shor-C0007-N12
- CERT-L0037-shell-diff-C0007-N24
- CERT-L0038-band-svd-C0007-N24
- CERT-L0039-sparse-support-C0007-N24
- CERT-L0040-sym-shor-C0007-N24
- CERT-L0041-shell6-sym-C0007-N24
- CERT-L0042-sym-strengthen-techo-C0007-N24
- CERT-L0043-onepol-shellblock-C0007-N24
- CERT-L0044-physical-fullsym-C0007-N24
- CERT-L0045-streaming-fullsym-C0007-N24
- CERT-L0046-greedy-fullsym-C0007-N24
- CERT-L0047-greedy-fullsym-C0007-N24
- Lean: `EnstrophyEnergy`, `EnstrophyEnergyRat`, `ComparisonODE`, `Certificates`

## Rejected

- C-0001, C-0002 (Stokes max-mode Ω(0.02)≈27.78 > M≈20.32), C-S-0001

## Explicit non-goals for the next sprints

- Re-deriving cancellation ODEs to “close” bounds already settled by L-0018
- Claiming continuum regularity / blow-up from Galerkin constants (`∼M/ν² → ∞`)
- Proposing all-IC M below the Stokes floor at the claimed T
