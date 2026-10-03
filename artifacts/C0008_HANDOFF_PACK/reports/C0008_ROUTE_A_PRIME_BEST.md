# C-0008 Route A' — Cluster + min(F,1-inf)

Mode **full**, 12 clusters, runtime 4166.68s.

| Quantity | Value |
|----------|-------|
| I_term^hi (Route A' best) | 37.9196 |
| Phase E per-shell I_hi | 170.2216 |
| Ratio A'/Phase E | 0.2228 |
| Omega(T)^hi | 54.2312 |
| Closes C-0008 | False |

Certificate: `certificates/CERT-L0073-C0008-route-a-prime-best.json`

Verify:
```bash
python verify_c0008_cluster_certificate.py certificates/CERT-L0073-C0008-route-a-prime-best.json experiments/terminal_weighted/shell_manifest_cluster_best_full.json
```

Per-cluster min(per-shell best sum, B_C max f(r)); B_C from best cluster combined bound. Triangle across clusters remains.
