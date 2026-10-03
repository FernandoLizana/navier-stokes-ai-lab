# C-0007 / C-0008 Structured Closure Map (L-0056)

> Finite Galerkin T^3, N<=24, dealias 2/3. **Not continuum. Not Clay.**

- **All-IC C-0007:** PROVED (M~41.7434, L-0070/L-0024)
- **Sharp C-0008:** OPEN (target M~41.2840)
- **Low-slab stretch target C_dagger:** 9.562010
- **Proved restricted subclasses:** 0 (C-R-0002..0014)

## Stretch constant ladder (attack C_dagger)

| Tier | Subclass / method | C bound | vs C_dagger | Evidence |
|------|-------------------|---------|-------------|----------|
| lower | rank-1 full dealias | ~3.45 | yes | N7 lower L-0049 |
| **best N5 SOS** | Hermitian {1,2,3} | **0.0000** | yes | CERT-L0052 |
| N5 SOS | Hermitian {1..4} | 0.894 | yes | CERT-L0054 |
| N5 SOS | Hermitian {1..5} | 1.305 | yes | CERT-L0058 |
| N5 SOS | Hermitian {1..6} | 1.701 | yes (interval) | CERT-L0060 |
| N5 SOS | Hermitian {1..8} | 1.852 | yes (interval COSMO) | CERT-L0061 |
| N5 SOS | Hermitian {1..9} | 2.223 | yes (interval COSMO) | CERT-L0062 |
| N5 SOS | Hermitian {1..10} | 2.440 | yes (interval COSMO) | CERT-L0062 |
| N5 SOS | Hermitian {1..12} | 2.783 | yes (interval COSMO) | CERT-L0062 |
| N5 SOS | one-pol C-R-0014 (12 shells) | 4.076 | yes | CERT-L0067 |
| N5 SOS | one-pol C-R-0013 (11 shells) | 4.015 | yes | CERT-L0065 |
| N5 SOS | one-pol C-R-0008 | 3.356 | yes | CERT-L0053 |
| Shor | greedy 41 shells C-R-0012 | 9.504 | yes | CERT-L0047 |
| extrap SOS | greedy 41 (cluster) | 5.39 | yes (extrap) | L-0055/L-0057 |
| techo Shor | all-IC fullsym | 25.93 | **no** | L-0048 |
| techo SOS | brute band extrap | ~14.7 | **no** | L-0051 |

## Proved restricted map (C-R)

| ID | C_fullsym | closes C_dagger | closes C-0007 subclass |
|----|-----------|-----------------|------------------------|

## Viable paths (honest)

1. **Done on laptop:** N5 SOS one-pol C-R-0014 (L-0064/67/68) + bands through shell 12.
2. **Proved without SOS:** C-R-0002..0014 via Shor/streaming (largest Hermitian: C-R-0012).
3. **Deferred:** greedy D=1772 SOS cluster (optional; C-R-0012 already Shor-closed).
4. **Abandoned:** brute Hermitian band SOS to full dealias (extrap ~14.7 > C_dagger).
5. **All-IC C-0007:** closed at L-0024 majorant (L-0070). **C-0008 sharp** still open.
6. **Refuted (2026-07-31):** L-0072 terminal shell integral (Phase D/E full-dealias); I_hi~78.8 vs I_*~1.30. See `conjectures/active/L-0072.json`.

## Terminal route (L-0071 / L-0072)

| Method | Route D I_hi | Closes C-0008? |
|--------|--------------|----------------|
| Frobenius Phase D | ~112.8 | No |
| Best (Phase E) | ~78.8 | No |
| Verifier | FAIL | Route refuted |

## N5 certificates (SOS Gram)

- CERT-L0052-sos-gram-band123-N24
- CERT-L0053-sos-gram-onepol-cr0008-N24
- CERT-L0054-sos-gram-band1234-N24
- CERT-L0058-sos-gram-band12345-N24
- CERT-L0060-sos-gram-band123456-N24
- CERT-L0061-sos-band12345678-N24
- CERT-L0062-sos-band123456789-N24
- CERT-L0062-sos-band1to10-N24
- CERT-L0062-sos-band1to12-N24
- CERT-L0065-sos-onepol-greedy-second-N24
- CERT-L0067-sos-onepol-greedy-second-r24-N24

## Lean N8 (L-0059)

- NSGalerkin.SosGram: native_decide C_ub_hi < C_dagger for N5 SOS certs
