# L-0073R: Cluster residual reweight (repair lemma)

> **Status:** candidate — replaces invalid `max(phi)*B_C` cluster bound  
> **Parent:** L-0071  
> **Historical L-0073:** `invalid_pending_cluster_reweight_proof` — I_hi≈37.92 is **not** rigorous

## Invalid bound (retired)

The Route A' cluster certificate used:

\[
\left\|\sum_{r\in C} \phi_r(t)\, M_r\right\|_{\mathrm{op}}
\;\stackrel{?}{\le}\;
\max_{r\in C}\phi_r(t)\cdot B_C
\]

where \(B_C\) bounds \(\|\sum_{r\in C} M_r\|\). This is **not valid** for general matrices:
cancellation inside \(\sum_r M_r\) is not captured by scalar \(\max\phi\).

## Valid identity

For any scalar \(\alpha_C\):

\[
\sum_r \phi_r M_r
= \alpha_C \sum_r M_r + \sum_r (\phi_r - \alpha_C) M_r.
\]

Triangle inequality on operator norm gives:

\[
\left\|\sum_r \phi_r M_r\right\|
\le
|\alpha_C|\,\left\|\sum_r M_r\right\|
+ \sum_r |\phi_r - \alpha_C|\,\|M_r\|.
\]

Certified upper bound (MPFR upward):

\[
\| \sum_r \phi_r M_r \|_{\mathrm{op}}
\le
|\alpha_C|\, B_C + \sum_{r\in C} |\phi_r - \alpha_C|\, B_r
\]

where \(B_r\) certifies \(\|M_r\|\) (one-inf baseline during repair) and
\(B_C\) certifies \(\|\sum_{r\in C} M_r\|\) via **residual combined pass** (future).

## Alpha candidates

- central shell weight \(\phi_{r_{\mathrm{center}}}\)
- \(\min_r \phi_r\), \(\max_r \phi_r\)
- midpoint of sorted \(\phi_r\)
- rigorous 1D optimizer over \(\alpha\in[\min,\max]\) (deferred)

## Evidence level

N4 until combined \(B_C\) uses exact coefficient enclosure and one-inf only baseline.

**Audit (2026-08-06):** `reports/L0073R_CLUSTER_AUDIT.md` — implementation bug (extra `max_phi` factor) fixed; corrected I_hi ≈ **38.52** ≈ historical ~37.92 numerically but now matches lemma identity. Still **N4 sketch**; triangle across clusters not N7.

## Does NOT close C-0008

Even with this lemma, sharp \(M\approx 41.28\) remains open.
