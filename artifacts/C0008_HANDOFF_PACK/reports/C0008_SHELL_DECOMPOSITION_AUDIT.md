# C-0008 — Shell decomposition audit (Phase D)

## Convention

Terminal-weighted operator at time `t`:

\[
\mathcal M(t) = \sum_{r \in \mathcal R} e^{-2\nu r (T-t)} \mathcal M_r,
\]

where each block `\mathcal M_r` is built with:

- output restricted to dealias shell `r`;
- inner product weight `w_r(t) = r e^{-2\nu r (T-t)}` at the output mode (L-0071);
- same Hermitian basis and 2/3 dealias triads as L-0048 fullsym.

At `t = T`: `w_r(T) = r`, matching the L-0048 terminal scan.

## Reconstruction (band N2)

On radii `{1..6}` at `n=24`:

- streaming `M(T)` Frobenius norm matches direct fullsym within band validation;
- per-shell blocks sum conservatively below full Frobenius (triangle on entries);
- `validate_direct_vs_tensor` passes on band (physical cubic bound, not Shor matricization).

Evidence level: **N2** for reconstruction; **N5** for certified Frobenius per shell.

## Full-dealias manifests (Phase D/E complete)

- `n = 24`, `D = 6748`, **87** non-empty dealias shells (`all_dealias_radii(24)`; IDs are radii, not 1..87).
- Band smoke test (6 shells, `D=160`): `experiments/terminal_weighted/shell_manifest.json`
- Phase D Frobenius (authoritative for Frobenius route): `experiments/terminal_weighted/shell_manifest_frobenius_full.json`
- Phase E best = min(Frobenius, 1-inf) per shell (authoritative for certificate): `experiments/terminal_weighted/shell_manifest_best_full.json`
- Certificate: `certificates/CERT-L0072-C0008-terminal-full-dealias.json`
- Reproduce: `.\reproduce_c0008_certificate.ps1` (verify-only) or `-Full` (recompute)

## L-0048 cross-check

| Quantity | Value | Role |
|----------|-------|------|
| L-0048 float `C(T)` | 25.925922… | N2 approximate `2√2 σ_max(M)` |
| Frobenius MPFR `C(T)` hi | see manifest | N5 rigorous **upper** bound |
| Classification | (a) float σ_max | **Not** certified upper |

## Limitations

1. `\|\mathcal M(t)\|_{op} \le \sum_r \phi_r(t) \|\mathcal M_r\|_{op}` — shell triangle may be loose.
2. Frobenius `\|\mathcal M_r\|_F` further loosens each block.
3. Band reconstruction does **not** substitute for full-dealias all-IC proof.
