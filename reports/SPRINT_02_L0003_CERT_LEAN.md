# Sprint 02 — L-0003 N5 Certificate (k≤3) + Lean Skeleton

**Route:** B. **No Clay claim.**

## Priority choice

After C-S-0002 survived stress, the highest-value next step was to **raise
evidence**, not re-run the same float attack. Two complementary layers:

1. **N5:** certify the L-0003 uniform constant for `|k|≤3` (L-0006 algebraic
   fails there).
2. **Toward N8:** Lean 4 formalization of the elementary discrete identity
   `Ω ≤ K² E`.

## 1. Why not L-0006 algebraic for k≤3

| shell | `K²` | `M` | `a√Ω0·t` (hi) | closes? |
|---|---|---|---|---|
| `|k|≤2` | 4 | 33 | 0.398 | yes → CERT-L0006 |
| `|k|≤3` | 9 | 123 | **1.153** | **no** → fall back |

## 2. CERT-L0003-shell-k3 (N5)

| quantity | value |
|---|---|
| formula | `Ω(t) ≤ max(K²E0, 6ME0²/ν²)` for all `t≥0` |
| `Ω0_hi` | 4.5 |
| `Ω_eq_hi` | 18450 |
| **guaranteed bound** | **18450** |
| verifier | PASS |

Loose (`∼ M/ν²`) but **globally closing** under hard shell truncation.
Artifact: `certificates/CERT-L0003-shell-k3.json`.

## 3. Lean skeleton (`lean/NSGalerkin`)

Formalized and **built successfully** (`lake build`):

```
NSGalerkin.EnstrophyEnergy.enstrophy_le_K2_energy
  : boundedBy K2 xs → weightSum xs ≤ K2 * energySum xs
```

Pure Lean 4 core (no Mathlib). Discrete Nat-list identity only —
**not** continuum NS, **not** Clay. Evidence: N8 *skeleton* for one step.

## 4. Honesty audit

- N5 constants certify truncated-Galerkin formulas only.
- Lean proves a combinatorial inequality, not the PDE.
- C-S-0002 remains an unproved N6 conjecture about full-grid evolution.
- Nothing here addresses Clay Routes A–D.

## 5. Reproduce

```
make sprint2-l0003-cert
make lean-build
make test
```

## 6. Next highest-value candidates

- Certify L-0003 also for the **full** dealiased grid at N=12 (not just shell).
- Lean: add the comparison-ODE closure step for L-0003 (still finite Nat/Real).
- Avoid treadmill re-stress of C-S-0002 unless a new IC family appears.
