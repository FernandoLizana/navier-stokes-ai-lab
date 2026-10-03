# L-0073 — Route A' cluster-best integral bound

**Status:** **best known terminal bound** — C-0008 still **`exploring`**

**Parent:** L-0071 / L-0072 follow-up

**Evidence:** N5 cluster manifest + verifier FAIL (expected)

## Result (2026-07-31, full-dealias)

| Quantity | Route A' best (12 clusters) | Phase E per-shell |
|----------|----------------------------|-------------------|
| `I_term^hi` | **37.92** | 78.80 (Route D) / 170.2 (Route A) |
| `Ω(T)^hi` | **54.23** | 68.69 |
| Verifier | **FAIL** | FAIL |

**Does not close** M ≈ 41.284. Best terminal pipeline to date; gap ~29× on `I_*`.

## Method

Per cluster `C`:

`I_C ≤ min( Σ_{r∈C} B_r^best f(r),  B_C^best · max_{r∈C} f(r) )`

with `B_C^best = min(Frobenius, 1-inf)` on combined cluster output (single D² pass).

## Artifacts

- Manifest: `experiments/terminal_weighted/shell_manifest_cluster_best_full.json`
- Certificate: `certificates/CERT-L0073-C0008-route-a-prime-best.json`

## Verifier

```bash
python verify_c0008_cluster_certificate.py \
  certificates/CERT-L0073-C0008-route-a-prime-best.json \
  experiments/terminal_weighted/shell_manifest_cluster_best_full.json
```

**Expected:** FAIL (`I_term_hi >= I_star_lo`, `omega_T_hi >= M_target`).
