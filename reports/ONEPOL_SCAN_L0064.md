# One-pol structured scan (L-0064)

> Finite Galerkin T³, N≤24. **Not all-IC. Not Clay.**

- **C_dagger:** 9.562010
- **One-pol SOS ratio (L-0053):** 0.3959
- **Candidates:** 82 | **Closing:** 46

## Best closing witness

| Field | Value |
|-------|-------|
| Pattern | greedy one-pol growth |
| Branch | second |
| Shells | `[1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 19, 24]` |
| D | 178 |
| C_Shor | 9.313031 |
| Est SOS C_ub | 3.6872 |
| SOS tier | cosmo_laptop |

## Greedy one-pol growth

- **first:** `[1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 18]` (11 shells)
- **second:** `[1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 19, 24]` (12 shells)

## SOS laptop queue (new, D≤356)

| scan_id | D | C_Shor | est SOS | tier | radii |
|---------|---|--------|---------|------|-------|
| greedy-second | 178 | 9.3130 | 3.687 | cosmo_laptop | `1,2,3,4,5,6,8,10,12,16,19,24` |
| greedy-first | 166 | 9.5565 | 3.784 | cosmo_laptop | `1,2,3,4,5,6,8,10,12,16,18` |
| cr0008+s14-first | 140 | 8.9370 | 3.538 | clarabel_laptop | `1,2,3,4,5,6,8,14` |
| cr0008+s17-first | 140 | 8.6206 | 3.413 | clarabel_laptop | `1,2,3,4,5,6,8,17` |
| cr0008+s21-first | 140 | 8.4976 | 3.364 | clarabel_laptop | `1,2,3,4,5,6,8,21` |
| cr0008+s14-second | 140 | 9.1304 | 3.615 | clarabel_laptop | `1,2,3,4,5,6,8,14` |
| cr0008+s17-second | 140 | 8.9041 | 3.525 | clarabel_laptop | `1,2,3,4,5,6,8,17` |
| cr0008+s21-second | 140 | 8.8142 | 3.490 | clarabel_laptop | `1,2,3,4,5,6,8,21` |
| cr0008+s18-first | 128 | 8.5363 | 3.380 | clarabel_laptop | `1,2,3,4,5,6,8,18` |
| cr0008+s18-second | 128 | 8.7754 | 3.474 | clarabel_laptop | `1,2,3,4,5,6,8,18` |
| cr0008+s25-first | 122 | 8.4833 | 3.359 | clarabel_laptop | `1,2,3,4,5,6,8,25` |
| cr0008+s25-second | 122 | 8.7305 | 3.457 | clarabel_laptop | `1,2,3,4,5,6,8,25` |
| cr0008+s10-first | 116 | 9.0582 | 3.586 | clarabel_laptop | `1,2,3,4,5,6,8,10` |
| cr0008+s11-first | 116 | 8.8291 | 3.496 | clarabel_laptop | `1,2,3,4,5,6,8,11` |
| cr0008+s13-first | 116 | 9.4235 | 3.731 | clarabel_laptop | `1,2,3,4,5,6,8,13` |

## Consecutive ladder (first branch)

| k | D | C_Shor | closes | est SOS |
|---|-----|--------|--------|---------|
| 1 | 6 | 0.0000 | yes | 0.000 |
| 2 | 18 | 2.0000 | yes | 0.792 |
| 3 | 26 | 3.2404 | yes | 1.283 |
| 4 | 32 | 3.5355 | yes | 1.400 |
| 5 | 56 | 5.5197 | yes | 2.185 |
| 6 | 80 | 7.6843 | yes | 3.042 |
| 6 | 80 | 7.6843 | yes | 3.042 |
| 8 | 92 | 8.4764 | yes | 3.356 |
| 9 | 122 | 9.8490 | no | 3.899 |
| 10 | 146 | 10.4248 | no | 4.127 |
| 11 | 170 | 12.6327 | no | 5.002 |
| 12 | 178 | 12.7901 | no | 5.064 |
| 13 | 202 | 14.2396 | no | 5.638 |
| 14 | 250 | 16.1039 | no | 6.376 |
| 14 | 250 | 16.1039 | no | 6.376 |
| 16 | 256 | 16.1096 | no | 6.378 |
| 17 | 304 | 17.5176 | no | 6.936 |
| 18 | 340 | 19.2432 | no | 7.619 |
| 19 | 364 | 20.1134 | no | 7.963 |
| 20 | 388 | 20.8871 | no | 8.270 |
| 21 | 436 | 23.1053 | no | 9.148 |
| 22 | 460 | 23.3352 | no | 9.239 |
| 22 | 460 | 23.3352 | no | 9.239 |
| 24 | 484 | 23.8808 | no | 9.455 |
| 25 | 514 | 25.1500 | no | 9.957 |
