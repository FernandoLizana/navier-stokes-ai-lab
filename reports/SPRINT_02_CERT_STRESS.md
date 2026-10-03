# Sprint 02 — L-0006 Interval Certificate + C-S-0002 Stress

**Route:** B. **Claim status:** exploratory / finite-dimensional only.
**No Clay Millennium claim is made anywhere in this sprint.**

## 1. What was done

1. **Extended interval arithmetic** (`validation/intervals.py`): added
   outward-rounded division (`__truediv__`), square root (`sqrt`), and square
   (`sqr`) to the existing `Interval` type.
2. **First N5 certificate of the project** — `CERT-L0006-shell-k2`:
   a machine-checked interval enclosure of the L-0006 algebraic short-time
   bound for the hard shell `|k|≤2`, `N=12`, `E0=0.5`, `t∈[0,0.02]`.
3. **Independent verifier** (`l0006_certificate_verify.py`) that re-derives the
   exact shell data from first principles and rejects understated bounds.
4. **Bounded stress campaign against C-S-0002** (`cs0002_stress.py`):
   structured + random + CMA-double-polish starts, all behind integrator and
   resolution gates.
5. **Theorem dossier** consolidating L-0001…L-0006 with honest scope
   (`docs/THEOREM_DOSSIER_GALERKIN.md`).

## 2. Certificate CERT-L0006-shell-k2 (N5)

| quantity | value | source |
|---|---|---|
| `K²` (exact integer) | 4 | shell lattice on `N=12` |
| `M` (exact integer) | 33 | Euclidean ball radius 2 |
| `a√Ω0·t` (upper) | 0.39799… | interval, `<1` ⇒ closes |
| **guaranteed** `Ω(t) ≤` | **5.5186…** | outward enclosure |

Verifier result: **PASS** (all checks true). Tampering tests confirm it rejects
an understated bound and a wrong mode count.

> Interpretation: this rigorously certifies the *constant* of a bound for the
> **hard-truncated shell ODE**. It is inviscid-upper, so it holds for all `ν≥0`.
> It says nothing about the continuum PDE and nothing about Clay.

## 3. C-S-0002 stress result

| field | value |
|---|---|
| proposed `M` | 48.15928 |
| best **gated** `Ω` found | **34.41168** (`abc_n16`) |
| refuted? | **No** |
| runs | 12 (structured / random / CMA-double-polish) |
| campaign wall time | ≈ 210 s |

C-S-0002 **survived**. The best gated enstrophy rose from the prior support
`32.106` to `34.412` but remains well below `M≈48.16`. The active record was
updated (support refreshed; `M` **not** loosened).

## 4. Honesty audit

- The certificate is **N5** for a **finite** statement only. The underlying
  inequality remains an **N7** hand proof.
- C-S-0002 remains an **N6** conjecture about the **full-grid** evolution of
  shelled ICs; the L-000x lemmas bound only the **hard-truncated** system, so
  they do **not** prove it.
- No result here approaches, tests, or bears on the continuum regularity /
  blow-up dichotomy of Clay Routes A–D.

## 5. Reproduce

```
make sprint2-cert-stress      # builds+verifies certificate, runs C-S-0002 stress
make test                     # includes test_l0006_certificate.py
```

Artifacts: `certificates/CERT-L0006-shell-k2.json`,
`conjectures/active/C-S-0002.json`, `docs/THEOREM_DOSSIER_GALERKIN.md`.

## 6. Next candidate steps

- Push the N5 certificate to `|k|≤3` (check whether the algebraic branch still
  closes at `t=0.02`; if not, certify the L-0003 uniform ceiling instead).
- Begin a Lean 4 skeleton for the elementary step `Ω ≤ K²E` (finite sum), the
  cleanest candidate to move one statement from N7 toward N8.
- Continue bounded stress on C-S-0002 only if a qualitatively new IC family is
  introduced (avoid treadmill re-runs of the same battery).
