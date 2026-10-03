# Sprint 02 — N16 L-0003 certificate + Lean ℚ enstrophy lift

**Route:** B. **No Clay claim.**

## Order

1. **CERT-L0003-full-dealias-N16** (N5) — documents `∼M/ν²` growth.
2. **Lean `EnstrophyEnergyRat`** — structural `Ω ≤ K²E` over nonnegative rationals.

## 1. Growth of the uniform Galerkin cap

| N | K² | M | guaranteed Ω-cap | ratio vs N12 |
|---|---|---|---|---|
| 12 | 27 | 343 | 51450 | 1 |
| 16 | 75 | 1331 | **199650** | M and cap both ×≈3.88 |

Exact identity (Lean kernel + Python bridge): for `(E0,ν)=(½,1/10)`,
`Ω_eq = 150·M`. Cap ratio equals M ratio — the bound **diverges** as
resolution → ∞, as required by honesty about truncated Galerkin.

Artifact: `certificates/CERT-L0003-full-dealias-N16.json` (verifier PASS).

## 2. Lean ℚ lift (`EnstrophyEnergyRat`)

`enstrophy_le_K2_energy_rat`: for mode energies `eᵢ = nᵢ/D` (`D≠0`, `nᵢ∈ℕ`),
`Σ k²·e ≤ K²·Σ e` whenever every `k² ≤ K²`.

Proved by clearing `D` and reducing to the Nat lemma (no Mathlib). Matches
Fourier `|û_k|² ≥ 0` with a common denominator (w.l.o.g. for finite lists).

Lean also proves `full_N16_cap = 199650` and `full_growth_N12_to_N16`.

## 3. Reproduce

```
make sprint2-n16-rat
make lean-build
make test
```

## 4. Honesty

- N16 cert is still truncated Galerkin; growth `∼M/ν²` is a **feature of the
  estimate**, not continuum blow-up.
- Lean ℚ lemma is combinatorial / arithmetic, not PDE theory.
- No Clay Routes A–D implication.
