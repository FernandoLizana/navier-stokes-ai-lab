# Sprint — L-0010 Short-Time C-S Bound

**Route:** B. **No Clay claim.**

## Result

First positive theorem in the C-S line:

> **C-S-0002-SHORT (N7):** For N=16, IC `|k|_∞≤4`, energy 0.5,
> `Ω(t) ≤ 48.159…` for all `t ∈ [0, T*]` with **`T* ≈ 0.001011`**.

The original C-S-0002 horizon `T=0.02` remains **open** (`T*` is ~20× shorter).

## Method (L-0010)

Duhamel bootstrap: under `Ω≤R`, high-mode energy satisfies
`E_H ≤ T² (U_L + W√ε)² R` with `ε=(R-B)/K²`. Closes iff that is `≤ ε`.

Also recorded: inviscid two-scale `T_max ≈ 0.00359` (L-0009); viscosity does
not push the C-S regime to `T=0.02`.

## N5

`certificates/CERT-L0010-duhamel-CS0002-short-N16.json` certifies the bootstrap
at `T_cert = 0.99 T*`.

## Artifacts

- `conjectures/proved_restricted/L-0010.json`
- `conjectures/proved_restricted/C-S-0002-SHORT.json`
- `reports/SPRINT_02_L0010.md`

## Reproduce

```
make sprint2-l0010
make test
```
