# Sprint — Refute C-0002; prove C-0003 via L-0018

## Headline

**C-0002 is refuted (N2)** by a single dealias Stokes mode `|k|²=75`:

| Quantity | Value |
|----------|-------|
| C-0002 M | ≈ 20.320 |
| Ω_ETD(0.02) | **≈ 27.7807** |
| Ω Stokes exact | `75·0.5·e^{-0.3}` ≈ 27.7807 |
| ETD↔RK4 rel | ~10⁻¹¹ |
| Successor | **C-0003** with **M = 37.5**, **proved** by L-0018 |

Prior CMA stress (max≈14.27) only searched low-mode shells and missed this IC.

## C-0003 (proved_restricted, N7)

All-IC on N≤16, dealias Galerkin:

\[
\Omega(t)\le K^2(16)\,E_0 = 37.5\qquad\forall t\ge 0.
\]

Same lemma that closed C-S-0002 (shell IC is a subclass).

## Files

- `conjectures/rejected/C-0002.json`
- `conjectures/proved_restricted/C-0003.json`
- `certificates/CERT-L0018-envelope-C0003-N16.json`
- `datasets/candidates/c0002_stokes_k2_75_refuter.npz`
- `experiments/exploratory/sprint02_refute_c0002/summary.json`

## Reproduce

```bash
pytest python/ns_exploration/tests/test_c0002_refute.py python/ns_exploration/tests/test_l0018.py -q
make sprint2-refute-c0002
```

## Honesty

Finite Galerkin / dealias only. Not continuum. Not Clay.
No active all-IC or shell-IC enstrophy conjecture remains at these (N≤16, E0=0.5) parameters below the envelope.
