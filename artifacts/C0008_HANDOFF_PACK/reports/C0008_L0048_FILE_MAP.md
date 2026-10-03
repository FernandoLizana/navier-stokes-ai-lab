# L-0048 file map (terminal-weighted audit)

> Maps code that builds the physical Hermitian fullsym stretch tensor.  
> At \(t=T\), terminal-weighted construction must regress to this pipeline.

## Call chain for `C_fullsym_full = 25.925922025377965`

```
sprint02_l0048.main()
  → all_dealias_radii(24)
  → sparse_C_fullsym(24, radii)
      → build_sparse_fullsym_csr()
          → build_hermitian_basis(n, radii)     [sprintB_physical_tensor]
          → _iter_fullsym_updates()             [l0048_matrixfree_fullsym]
              → _inner_psi_N(b, m, Nloc)        [l0045_streaming_fullsym]
      → sparse_gram_opnorm_csr(M)
  → C = 2√2 · L_op
```

## Key files

| File | Functions |
|------|-----------|
| `python/ns_exploration/conjectures/l0048_fullsym_techo.py` | `lemma_l0048()`, `C_FULLSYM_FULL_DEALIAS` |
| `python/ns_exploration/conjectures/l0048_matrixfree_fullsym.py` | `sparse_C_fullsym`, `build_sparse_fullsym_csr`, `all_dealias_radii` |
| `python/ns_exploration/conjectures/l0045_streaming_fullsym.py` | `streaming_fullsym_M`, `_inner_psi_N`, `streaming_C_fullsym` |
| `python/ns_exploration/experiments/sprintB_physical_tensor.py` | `build_hermitian_basis`, `build_physical_G`, `z_to_uhat` |
| `python/ns_exploration/conjectures/l0040_sym_shor.py` | `_sym_col_index` |
| `python/ns_exploration/conjectures/l0021_quartic_shell.py` | `_leray_vec`, `shell_modes_by_r` |
| `python/ns_exploration/conjectures/l0023_cubic_dissipation.py` | `stretch_inner` (FFT cross-check) |
| `python/ns_exploration/spectral/fourier_conventions.py` | FFT normalization, Parseval |
| `python/ns_exploration/spectral/dealias.py` | `dealias_mask` |

## Terminal-weighted adapter (new)

| File | Role |
|------|------|
| `python/ns_exploration/terminal_weighted/weights.py` | `w_r(t)=r e^{-2νr(T-t)}` |
| `python/ns_exploration/terminal_weighted/tensor.py` | `_inner_psi_N_weighted`, `streaming_terminal_fullsym_M` |
| `python/ns_exploration/terminal_weighted/scan.py` | N2 temporal grid |

## Regression requirement

At `t=T`, `w_r(T)=r` ⇒ `_inner_psi_N_weighted` = `_inner_psi_N` ⇒  
`terminal_fullsym_bound(n, radii, T)` must match `streaming_C_fullsym(n, radii)`.

Full dealias sparse at `t=T` must match `25.925922025377965` (L-0048 frozen).

## Normalization checklist

- [x] FFT: `û = fftn(u)/N³`
- [x] Energy: `E = ½ Σ |û_k|²`
- [x] Enstrophy: `Ω = ½ Σ |k|² |û_k|²`
- [x] Hermitian: half-space + `(c,s)` per pol
- [x] Leray: `-P_s(Nloc)` after triad sum
- [x] Stretch weight at T: `|k|²` on output mode
- [x] full_sym factor `1/6` perms + global `×0.5` on M

*Finite Galerkin only.*
