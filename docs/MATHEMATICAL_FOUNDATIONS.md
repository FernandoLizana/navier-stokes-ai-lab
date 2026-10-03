# Mathematical Foundations

**Route focus:** B (periodic, \(f=0\)); secondary D.  
**Evidence:** N7 for classical identities cited below; no new theorems claimed here.

## 1. Incompressible Navier–Stokes on \(\mathbb{T}^3\)

Let \(\mathbb{T}^3 = \mathbb{R}^3/(2\pi\mathbb{Z})^3\) (NS-MRL computational convention; equivalent to \(\mathbb{R}^3/\mathbb{Z}^3\) after rescaling). Viscosity \(\nu>0\).

\[
\partial_t u + (u\cdot\nabla)u = -\nabla p + \nu\Delta u + f,
\qquad
\nabla\cdot u = 0,
\qquad
u|_{t=0}=u_0,\quad \nabla\cdot u_0=0.
\]

**Terms (never dropped without justification):**

| Term | Role |
|------|------|
| \(\partial_t u\) | Unsteady acceleration |
| \((u\cdot\nabla)u\) | Advection / inertial nonlinearity |
| \(-\nabla p\) | Pressure (enforces incompressibility) |
| \(\nu\Delta u\) | Viscous diffusion / dissipation |
| \(f\) | External force (zero on Route B) |

## 2. Leray projection

On \(\mathbb{T}^3\), the Helmholtz–Leray projector \(P\) maps to divergence-free fields. In Fourier modes \(k\in\mathbb{Z}^3\setminus\{0\}\):

\[
P_k = I - \frac{kk^\top}{|k|^2},\qquad
\widehat{Pu}_k = P_k\hat u_k.
\]

The projected form eliminates pressure:

\[
\partial_t u + P\bigl((u\cdot\nabla)u\bigr) = \nu\Delta u + Pf.
\]

**Zero mode:** \(k=0\) carries mean velocity; pressure has no effect on the mean if \(\widehat{\nabla p}_0=0\). Mean momentum evolves by mean force (Route B: conserved mean if \(f=0\)).

## 3. Vorticity formulation

\(\omega=\nabla\times u\):

\[
\partial_t\omega + (u\cdot\nabla)\omega = (\omega\cdot\nabla)u + \nu\Delta\omega + \nabla\times f.
\]

| Term | Name |
|------|------|
| \((u\cdot\nabla)\omega\) | Vorticity transport |
| \((\omega\cdot\nabla)u\) | **Vortex stretching** (absent in 2D) |
| \(\nu\Delta\omega\) | Viscous diffusion of vorticity |
| \(\nabla\times f\) | Forcing curl |

Vortex stretching is the principal obstruction to closing 3D enstrophy estimates.

## 4. Energy identity (smooth solutions, \(f=0\))

For smooth divergence-free \(u\) on \(\mathbb{T}^3\),

\[
\frac{d}{dt}\frac12\|u\|_{L^2}^2 + \nu\|\nabla u\|_{L^2}^2 = 0.
\]

The nonlinear term cancels after integration by parts and \(\nabla\cdot u=0\) (Leray–Hopf energy equality for smooth solutions). Enstrophy \(\mathcal{E}=\frac12\|\omega\|_{L^2}^2\) is **not** a priori controlled in 3D.

## 5. Scale invariance (reminder)

The NS scaling \(u_\lambda(x,t)=\lambda u(\lambda x,\lambda^2 t)\) leaves the equation invariant (\(f=0\)). Critical spaces (e.g. \(\dot H^{1/2}\), \(BMO^{-1}\)) are those invariant under this scaling. Clay Routes A/B ask for **all** smooth data, which may be large in critical norms.

## 6. Local existence (classical, cited)

For smooth divergence-free data on \(\mathbb{T}^3\) or \(\mathbb{R}^3\), local-in-time unique smooth solutions exist; the maximal time is characterized by blow-up criteria (Beale–Kato–Majda and refinements). See `docs/REGULARITY_CRITERIA.md` and `literature/`.

## 7. What this document does not claim

- Global regularity or blow-up.
- Any new estimate.
- That numerical spectral decay implies smoothness of the continuum solution.

## Primary references

- J. Leray, *Acta Math.* (1934).
- E. Hopf, *Math. Nachr.* (1951).
- C. Fefferman, Clay official description (2000/2006).
- P. Constantin & C. Foias, *Navier–Stokes Equations*, Chicago (1988).
- R. Temam, *Navier–Stokes Equations*, North-Holland (1984).
