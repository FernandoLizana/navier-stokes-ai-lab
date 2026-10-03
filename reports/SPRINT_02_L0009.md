# Sprint — L-0009 Two-Scale Embedding

**Route:** B. **No Clay claim.**

## Idea

Split `||u||_∞` into shell + high-mode pieces (`M_L`, `M_H`) and solve the
comparison ODE for high-mode energy starting at zero:
`z' ≤ (U_L + W z) √(B + K_full² z²)`, closed form via sinh substitution.

## Results

| case | L-0009 | L-0008 | L-0007 |
|---|---|---|---|
| N=12, eucl. ≤2 | **103.58** | 394 | 1573 |
| N=16, ℓ∞≤4 (C-S) | no close | 199650 | 199650 |

Improvement vs L-0008 at the N=12 model problem: **~3.8×**.
C-S-0002 still open (two-scale integral saturates before T=0.02 on N=16).

## N5

`CERT-L0009-twoscale-euclidean-k2-N12`: Ω ≤ **103.584** (verifier PASS,
including `F(u_hi)-F(0) ≥ K t`).

## Reproduce

```
make sprint2-l0009
make test
```

## Bottleneck now

On the C-S-0002 regime the two-scale comparison ODE ceases to exist past a
finite time < 0.02. Need either a shorter horizon, a viscosity-aware estimate,
or a sharper high-mode production bound (not just `√(3 M_H) √(2 E_H)`).
