# L-0071: Terminal-weighted enstrophy identity (candidate)

> **Status:** candidate lemma, N7 derivation + N2 engineering tests.  
> **Does NOT close C-0008.** Finite Galerkin T³, N≤24, dealias 2/3 only.

## 1. Setup

Truncated NS Galerkin:
\[
\dot u = -\nu A u - N(u), \qquad A = -\Delta \text{ on retained modes}.
\]

Energy \(E(t)=\tfrac12\|u\|_2^2\), enstrophy \(\Omega(t)=\tfrac12\langle u,Au\rangle\).

Parameters frozen: \(\nu=0.1\), \(E_0=0.5\), \(T=0.02\), \(N\le24\).

## 2. Weight operator

\[
W(t) = A e^{-2\nu A(T-t)}.
\]

On each shell \(r=|k|^2\): eigenvalue \(w_r(t) = r e^{-2\nu r(T-t)}\).

## 3. Terminal functional

\[
\Phi(t,u) = \tfrac12 \langle u, W(t) u \rangle
= \tfrac12 \sum_k |k|^2 e^{-2\nu|k|^2(T-t)} |\hat u_k|^2.
\]

## 4. Derivative (viscous cancellation)

Since \(W'(t) = 2\nu A W(t)\),
\[
\frac{d}{dt}\Phi(t,u(t))
= \langle W(t)u, \dot u \rangle + \tfrac12 \langle u, W'(t) u \rangle
\]
\[
= \langle W(t)u, -\nu Au - N(u) \rangle + \nu \langle u, A W(t) u \rangle
= -\langle W(t)u, N(u) \rangle.
\]

## 5. Terminal and initial values

At \(t=T\): \(W(T)=A\), hence \(\Phi(T,u(T))=\Omega(T)\).

At \(t=0\), using \(E(0)\le E_0\):
\[
\Phi(0,u_0) \le E_0 \max_{r\in\mathcal R_N} r e^{-2\nu r T}.
\]
For N=24 dealias this bound equals the Stokes floor `40.82462307930434` (certified hi in code).

## 6. Integral inequality

Integrating from 0 to T:
\[
\Omega(T) = \Phi(0,u_0) - \int_0^T \langle W(t)u, N(u)\rangle\,dt
\le \text{Stokes floor} + \int_0^T |\langle W(t)u, N(u)\rangle|\,dt.
\]

Using \(\|u(t)\|_2\le 1\) and homogeneity:
\[
|\langle W(t)u, N(u)\rangle| \le \frac{c_{\mathrm{term}}(t)}{2\sqrt2},
\]
where
\[
c_{\mathrm{term}}(t) = 2\sqrt2 \sup_{\|z\|=1} |\langle W(t)z, N(z)\rangle|.
\]

Hence
\[
\Omega(T) \le \text{Stokes floor} + \frac{1}{2\sqrt2}\int_0^T c_{\mathrm{term}}^{ub}(t)\,dt.
\]

## 7. C-0008 threshold (not proved)

Sufficient condition for target \(M=41.283999509509286\):
\[
\int_0^T c_{\mathrm{term}}^{ub}(t)\,dt \le 1.2993127556607498.
\]

## 8. Normalizations (audit)

- Fourier: \(\hat u_k = \mathrm{fftn}(u)/N^3\), Parseval \(E=\tfrac12\sum|\hat u_k|^2\).
- Hermitian real fields: half-space \(k\) with \((c,s)\) coords per Leray pol (see `sprintB_physical_tensor`).
- Nonlinearity: triad sum with factor `1j * (va·kq)`, then Leray \(-P_s(\cdot)\).
- Stretch observable at \(t=T\) matches L-0048/L-0045 `_inner_psi_N` with weight \(|k|^2\).

## 9. Evidence level

| Claim | Level |
|-------|-------|
| Identity \(\Phi'=-\langle Wu,Nu\rangle\) | N7 symbolic |
| Stokes floor bound | N7 + N5 hi in code |
| \(c_{\mathrm{term}}^{ub}\) via fullsym | N2 until interval cover |
| C-0008 closure | **not claimed** |

*Not continuum. Not Clay.*
