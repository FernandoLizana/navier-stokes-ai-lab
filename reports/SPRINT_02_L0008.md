# Sprint — L-0008 Cascade Bound (`E_H(0)=0`)

**Route:** B. **No Clay claim.**

## Idea

L-0007 only used shell `Ω0`. L-0008 uses that high-mode energy starts at zero:

`Ω ≤ B + K_full² E_H`, `E_H(0)=0`, and a comparison ODE for `z=√E_H` yields

`Ω(t) ≤ B cosh²(A K_full t)`, `A = √2 √(3M_full) √E0`, `B = K_IC² E0`.

## Results

| case | L-0008 | L-0007 | L-0003 |
|---|---|---|---|
| N=12, eucl. ≤2 | **394** | 1573 | 51450 |
| N=12, ℓ∞ ≤4 | 4730 | 18870 | 51450 |
| N=16, ℓ∞ ≤4 (C-S) | explodes | 199650 | **199650** |

- ~4× sharper than L-0007 at N=12 short time.
- C-S-0002 class still unproved (gap ~4146× to M≈48).

## N5 certificate

`CERT-L0008-cosh-euclidean-k2-N12` — verifier PASS.

## Reproduce

```
make sprint2-l0008
make test
```

## Honest next bottleneck

The factor `√(3 M_full)` treats every retained mode as if it could sit in the
`||u||_∞` sum at full strength. Closing C-S-0002 needs a sharper embedding or
a Gevrey/analyticity radius that keeps effective `M` small on `[0,0.02]`.
