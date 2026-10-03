# NSGalerkin — Lean 4 skeleton (finite-mode identities)

**Purpose.** Formalize elementary combinatorial steps used by the L-0001…L-0006
Galerkin bounds, starting with the discrete enstrophy–energy inequality.

**Honesty.** This package proves statements about finite lists of natural
numbers. It does **not** formalize continuum Navier–Stokes, and it has **no
Clay Millennium implication**. Evidence target: N8 *skeleton* for one lemma.

## Build

```bash
cd lean/NSGalerkin
lake build
```

## Contents

| module | statement |
|---|---|
| `NSGalerkin.EnstrophyEnergy` | `Σ k²·e ≤ K² · Σ e` when every `k² ≤ K2` (Nat) |
| `NSGalerkin.EnstrophyEnergyRat` | same over ℚ≥0 with common den `D` (`eᵢ=nᵢ/D`) |
| `NSGalerkin.ComparisonODE` | sign of L-0003 comparison force; `max(Ω0,Ω_eq)` cap; `Ω_eq·ν² = 6ME0²` |
| `NSGalerkin.Certificates` | exact CERT-L0003 caps (shell-k3, N12, N16) + growth witness |

## Relation to Python certificates

- Python N5 certificates enclose the *numeric constants* of L-0003 / L-0006
  (`certificates/CERT-*.json`), including full-dealias N=12 and N=16.
- This Lean package aims at N8 for the *structural* inequalities behind those
  bounds, plus exact ℚ arithmetic of the certificate constants.
- Together they cover different layers of the same truncated-Galerkin stack.
