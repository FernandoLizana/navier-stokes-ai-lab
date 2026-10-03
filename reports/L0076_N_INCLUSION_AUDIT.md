# L-0076 N≤24 Inclusion Audit

**Evidence level:** N7_sketch_empirical  
**C-0008:** exploring

## Wavevector embedding (n=24 majorant)

- All subsets: `True`

## Band majorant ladders

| Band | Monotone (n=12→24) |
|------|---------------------|
| {1,2,3} | True |
| {1..6} | True |

## Repair full-dealias n=24 (Route A)

| Field | Value |
|-------|-------|
| Shells | 87 |
| Complete | True |
| I_hi | 170.5871515665858369747011687095425020673 |

## Conditional conclusion

IF mode embedding + majorant monotonicity hold, then n=24 repair full-dealias I_hi majorizes embedded n'<=24 truncations.

## Formal gaps (still open)

- Embedding z_{n'} into z_n preserving terminal functional (N7).
- Per-n full-dealias certificates for all n<=24 (only n=24 complete).
- Monotonicity of sup f_n under mode-set inclusion not proved.
