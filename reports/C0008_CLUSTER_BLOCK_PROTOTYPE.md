# C-0008 — Cluster block bounds prototype (band)

> Prototype after L-0072 per-shell route refutation. **Band only (N2/N5 pilot).**

## Idea

Route A uses triangle inequality on shells:

\[
\|M(t)\|_{op} \le \sum_r \phi_r(t)\,\|M_r\|_{op}.
\]

For a cluster \(C = \{r_1,\ldots,r_k\}\), bound \(\|\sum_{r\in C} M_r\|_F\) directly
instead of \(\sum_{r\in C}\|M_r\|_F\). By Frobenius triangle,

\[
\Big\|\sum_{r\in C} M_r\Big\|_F \le \sum_{r\in C}\|M_r\|_F,
\]

so the cluster norm is **never worse** at \(t=T\).

## API

```python
from ns_exploration.terminal_weighted.opnorm_certified import (
    certified_cluster_frobenius_hi,
    compare_shell_vs_cluster_band,
)

cmp = compare_shell_vs_cluster_band(24, tuple(range(1, 7)), prec=128)
# cmp["cluster_tighter_ratio"] < 1 ⇒ cluster beats per-shell sum at t=T
```

## Band pilot result (n=24, shells 1..6)

| Quantity | Value |
|----------|-------|
| Sum per-shell `C_term^hi` | ~121.56 |
| Cluster `{1..6}` `C_term^hi` | ~56.80 |
| **Ratio (cluster/sum)** | **~0.47** |

Cluster Frobenius cuts the band t=T bound by ~53% vs triangle sum. Integral Route A
still uses per-shell temporal weights; a cluster-aware certificate is required to
transfer this gain to `I_term^hi`.


1. Integral Route A uses **per-shell** weights \((1-e^{-2\nu r T})/(2\nu r)\); cluster
   bounds need a reweighted temporal certificate (not implemented).
2. Full dealias cluster pass at D=6748 is deferred (cost similar to per-shell Frobenius).
3. Does **not** by itself close C-0008; Phase E best I_hi~78.8 still FAIL.

## Next steps

- Cluster-aware Route A' with per-cluster temporal weights
- Sparse Gram eig on terminal tensor (`tensor.py::terminal_sparse_at_T`)
- Hybrid: cluster inner shells, single shells on tail radii
