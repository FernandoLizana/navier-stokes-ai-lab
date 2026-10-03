# Sprint — L-0012 NS cancellations → longer C-S-0002-SHORT

## Result

| bound | T* for R≈48.16 (N=16, ℓ∞≤4) |
|---|---|
| L-0010 (Duhamel) | ≈0.001011 |
| L-0011 (div-free + viscous τ) | ≈0.001240 (~1.23×) |
| **L-0012** (NS cancellations for E_H) | **≈0.002284** (~1.84× vs L-0011) |

C-S-0002-SHORT updated to this horizon. N5: `CERT-L0012-cancel-CS0002-short-N16`.

Parent **C-S-0002** at T=0.02 remains **open** (L-0012 majorant at 0.02 ≫ R).

## Math (sketch)

High-mode energy after Leray + div-free cancellations:
```
dE_H/dt ≤ -2νκ_H E_H + 2 U_L √B √E_H + 2 W √B E_H
```
Linear comparison for `z=√E_H`, `z(0)=0` → closed Ω bound on short time.

## Artifacts

- `python/ns_exploration/conjectures/l0012_cancellation.py`
- `python/ns_exploration/validation/l0012_certificate.py`
- `conjectures/proved_restricted/L-0012.json`
- `conjectures/proved_restricted/C-S-0002-SHORT.json`
- `certificates/CERT-L0012-cancel-CS0002-short-N16.json`

## Commands

```
pytest python/ns_exploration/tests/test_l0012.py -q
make sprint2-l0012
```

## Honesty

Finite Galerkin, short horizon only. Not continuum. Not a Clay claim.
