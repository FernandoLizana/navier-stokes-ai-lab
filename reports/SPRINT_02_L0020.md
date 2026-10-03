# Sprint — L-0020 Duhamel H¹ (N_crit for C-0005)

## Result

Stokes–Duhamel in Ḣ¹ on N=24 dealias:

\[
\Omega(T)\le \tfrac12\bigl(\sqrt{2 S(T)E_0}+N_* I_\sigma\bigr)^2.
\]

| quantity | value |
|----------|-------|
| Stokes floor `S E0` | ≈40.825 |
| `I_σ=∫σ` | ≈0.210 |
| Young `N_*` | ≈78.69 (`ρ_★=6192`) |
| Young ⇒ `Ω_H1` | ≈327 (**worse than envelope 73.5**) |
| **N_crit for C-0005** | **≈9.666** |
| Empirical `‖N‖` | ≲2 |

**C-0005 not closed.** Path forward: certify any `N_* ≤ N_crit` (room vs empirics ≈4.8×).

## Files

- `conjectures/proved_restricted/L-0020.json`
- `certificates/CERT-L0020-duhamel-H1-N24.json`
- updated `conjectures/active/C-0005.json` notes

## Reproduce

```bash
pytest python/ns_exploration/tests/test_l0020.py -q
make sprint2-l0020
```

## Honesty

Finite Galerkin only. Young instantiation does not beat L-0018. Not Clay.
