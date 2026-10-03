# Sprint — Hygiene + L-0007 (shell IC → full grid)

**Route:** B. **No Clay claim.**

## 1. Hygiene

- Removed refuted records from `conjectures/active/` (`C-0001`, `C-S-0001`);
  they remain under `conjectures/rejected/`.
- Active set is now only **C-0002** and **C-S-0002**.
- `docs/ROADMAP.md` rewritten to match actual status (adversarial + N5/Lean
  no longer “deferred/stubs”).

## 2. L-0007 — the missing bridge

Shell-IC evolution on the **full** dealias grid (C-S class):
use full-grid stretch `(K_full, M_full)` but shell initial ceiling
`Ω0 ≤ K_IC² E0`.

| N | IC class | best proved | method | proves C-S-0002? |
|---|---|---|---|---|
| 12 | eucl. ≤2 | **1573** | L-0007-exp | no |
| 12 | ℓ^∞ ≤4 | 1.89e4 | L-0007-exp | no |
| 16 | ℓ^∞ ≤4 | 1.997e5 | L-0003-full | **no** (gap ≈4146×) |

So: L-0007 is a real improvement at N=12 short time, and at the C-S-0002
setting (N≤16, ℓ^∞≤4) it **quantifies** that current a priori technology is
~4000× too weak for M≈48.

## 3. N5 certificate

`CERT-L0007-exp-euclidean-k2-N12`: outward enclosure
`Ω(t) ≤ 1572.53` for `t∈[0,0.02]`. Verifier PASS.

## 4. Reproduce

```
make sprint2-l0007
make test
```

## 5. What this does *not* do

- Does not prove or refute C-S-0002.
- Does not approach continuum Clay Routes A–D.
- Next mathematical need: a sharper estimate that exploits that high modes
  start at *zero* (not just shell Ω0), e.g. bootstrap / Gevrey / short-time
  cascade control — still finite-dimensional first.
