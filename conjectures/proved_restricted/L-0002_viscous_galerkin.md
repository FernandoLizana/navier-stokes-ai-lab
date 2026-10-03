# L-0002 — Viscous Galerkin enstrophy comparison ODE

**Status:** proved_restricted (finite Fourier–Galerkin, mean-zero)  
**Evidence:** N7  
**Clay:** none

## Statement

With \(C = 2\sqrt{2}\sqrt{3M}\), \(a = C/2\), \(\Omega(0)\le K^2 E_0\), \(z_0=\sqrt{\Omega(0)}\),

\[
\frac{d\Omega}{dt} \le -2\nu\Omega + C\Omega^{3/2}
\]

via \(\|\nabla u\|_\infty \le \sqrt{3M}\sqrt{2\Omega}\) and \(\|\nabla\omega\|^2 \ge 2\Omega\).
The comparison ODE for \(z=\sqrt{\Omega}\) has solution

\[
z(t) = \frac{\nu}{a - (a - \nu/z_0)e^{\nu t}}
\]

while the denominator stays positive; else the *estimate* fails to close (not a proof of singularity).

## Caps on NS-MRL grids (E0=0.5, t=0.02, ν=0.1)

| N | L-0001 cap | L-0002 cap / status |
|---|------------|---------------------|
| 8 | 8.779820e+01 | ESTIMATE FAILS (t*=0.014918242015779914) |
| 12 | 1.061457e+04 | ESTIMATE FAILS (t*=0.006001236461296834) |
| 16 | 1.203975e+11 | ESTIMATE FAILS (t*=0.0018275084614762527) |

Auditor note: failure of L-0002 to close at large M is expected and must not be misread as Galerkin blow-up.
