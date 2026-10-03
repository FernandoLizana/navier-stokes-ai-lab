# Sprint 02 — Full-dealias L-0003 N5 + Lean comparison ODE

**Route:** B. **No Clay claim.**

## Order chosen

1. **First:** CERT-L0003 on the full 2/3 dealias grid (N=12) — extends the
   existing N5 pipeline immediately.
2. **Second:** Lean `ComparisonODE` — algebraic closure steps of the L-0003
   comparison force / uniform cap (N8 skeleton).

## 1. CERT-L0003-full-dealias-N12 (N5)

| quantity | value |
|---|---|
| support | full 2/3 dealias (`|k_i| < N/3`) |
| `K²`, `M` | **27**, **343** |
| `Ω0_hi` | 13.5 |
| `Ω_eq_hi` | 51450 |
| **bound** | **51450** (all `t≥0`) |
| verifier | PASS |

Artifact: `certificates/CERT-L0003-full-dealias-N12.json`.
Shell cert `CERT-L0003-shell-k3` regenerated with explicit `support` field.

## 2. Lean `ComparisonODE` (builds)

`lake build` OK. Formalized:

- `force_nonpos_of_dissip` / `force_nonneg_of_stretch` — sign of `a z²` vs `b z³`
- `equilibrium_of_exact_div` — exact Nat equilibrium
- `uniform_cap_closes` — `max(Ω0, Ω_eq)` + product form `Ω_eq·ν² = 6 M E0²`

Honesty: **Nat algebra only**; the differential inequality remains hand N7.
Does not formalize continuum NS.

## 3. Reproduce

```
make sprint2-full-l0003-lean
make lean-build
make test
```

## 4. Honesty audit

- Full-dealias bound is still truncated Galerkin; `Ω_eq ∼ M/ν² → ∞` as `N→∞`.
- Lean does not prove the PDE comparison principle — only the arithmetic that
  sits under CERT-L0003.
- No Clay Routes A–D implication.
