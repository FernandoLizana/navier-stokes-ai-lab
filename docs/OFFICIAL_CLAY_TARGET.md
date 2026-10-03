# Official Clay Mathematics Institute Target

**Project:** Navier–Stokes Millennium Research Laboratory (NS-MRL)  
**Evidence level of this document:** N7 (editorial synthesis of the official problem statement; not a new theorem)  
**Primary computational route:** **B**  
**Secondary route:** **D**  
**Routes never mixed in a single claim.**

This file records the four statements accepted by the Clay Mathematics Institute (CMI) official problem description, then explains them mathematically. The authoritative source is Fefferman’s official formulation for the Clay Millennium Problems (see literature matrix entry `Fefferman2000`).

---

## Literal routes (official formulation)

The official problem asks for a proof of **one** of the following.

### Route A — Global regularity on \(\mathbb{R}^3\), zero force

For smooth, divergence-free, rapidly decaying initial data \(u_0\) on \(\mathbb{R}^3\), and external force \(f \equiv 0\), there exists a smooth solution \(u(x,t)\) of the incompressible Navier–Stokes equations that exists for all \(t \geq 0\) and remains smooth.

### Route B — Global regularity on \(\mathbb{T}^3\), zero force

For smooth, divergence-free, **periodic** initial data \(u_0\) on the 3-torus \(\mathbb{T}^3 = \mathbb{R}^3/\mathbb{Z}^3\) (equivalently, \(2\pi\)-periodic after rescaling), and \(f \equiv 0\), there exists a smooth global-in-time solution.

### Route C — Blow-up on \(\mathbb{R}^3\) with smooth force

There exist smooth divergence-free initial data and a smooth external force on \(\mathbb{R}^3\) such that no smooth solution exists for all positive times (finite-time breakdown of smoothness).

### Route D — Blow-up on \(\mathbb{T}^3\) with smooth force

There exist smooth divergence-free **periodic** initial data and a smooth **periodic** force on \(\mathbb{T}^3\) such that no smooth global solution exists.

---

## Mathematical explanation of each route

### Common PDE

\[
\partial_t u + (u\cdot\nabla)u = -\nabla p + \nu\Delta u + f,\qquad
\nabla\cdot u = 0,\qquad
u(\cdot,0)=u_0,
\]
with \(\nu>0\). The pressure is determined (up to a constant) by the divergence-free constraint via the Leray projection.

### Route A vs B (regularity)

Both assert: **for every** admissible smooth \(u_0\) (rapidly decaying on \(\mathbb{R}^3\), or smooth periodic on \(\mathbb{T}^3\)), with \(f=0\), a unique smooth solution exists globally. The difficulty is **large data** and **possible uncontrolled vortex stretching** in 3D. Local existence for smooth data is classical; the Millennium question is **global** control.

### Route C vs D (breakdown)

Both assert existence of **some** smooth data (and smooth \(f\)) for which smoothness fails in finite time. On \(\mathbb{T}^3\), force is allowed in Route D; NS-MRL treats \(f\) as part of the constructed data and never hides singularities inside \(f\).

### What “smooth” and “physically reasonable” mean

In the official write-up, solutions must satisfy the PDE classically (or in a sense that implies classical smoothness after local existence theory), remain divergence-free, and obey the problem’s growth/decay or periodicity hypotheses. Weak Leray–Hopf solutions alone do **not** settle Routes A/B; non-uniqueness of weak solutions does **not** settle Routes C/D.

---

## NS-MRL primary choice: Route B

**Why \(\mathbb{T}^3\), \(f=0\):**

- Exact Fourier series \(u(x,t)=\sum_{k\in\mathbb{Z}^3}\hat u_k(t)\,e^{2\pi i k\cdot x}\) (or \(e^{ik\cdot x}\) on \([0,2\pi)^3\); **one** convention is fixed project-wide).
- Exact Leray projector per mode: \(P_k = I - kk^\top/|k|^2\) for \(k\neq 0\).
- Clean divergence-free constraint \(k\cdot\hat u_k=0\).
- Spectral tails, symmetries, and reproducible validated numerics.

**Secondary:** Route D (periodic blow-up with smooth \(f\)), in separate campaigns and certificates.

---

## What would constitute a solution

### For Route B (regularity)

A proof that for **all** smooth divergence-free periodic \(u_0\), \(\nu>0\), \(f=0\):

1. a solution exists for all \(t\geq 0\);
2. it remains smooth for all \(t\geq 0\);
3. it satisfies the exact PDE;
4. constants do not depend on an undemonstrated a priori bound;
5. the argument covers **large** data (not only small-data / perturbative regimes);
6. any Galerkin / approximation limit is justified.

### For Route D (breakdown)

An explicit smooth periodic divergence-free \(u_0\), explicit smooth periodic \(f\) (if used), a finite \(T>0\), and a rigorous proof that no smooth solution exists on \([0,T]\), with all errors controlled and any reduced profile reconstructed in the full 3D PDE.

---

## What would **not** constitute a solution

| Claim type | Why insufficient |
|------------|------------------|
| Floating-point blow-up of \(\|\omega\|_\infty\) | Mesh/aliasing/time-step artefacts (`≤ N3`) |
| Neural / PINN residual small in mean | No rigorous norm bound |
| Truncated Galerkin blow-up | Not the infinite-dimensional PDE |
| Regularity for one symmetry class | Not “for all data” on Route B |
| Small-data global regularity | Already classical in critical spaces; does not finish Clay |
| Weak non-uniqueness | Different statement |
| Setting \(\nu=0\) (Euler) | Different equation |
| Hyperviscosity / different dimension / walls | Outside official formulation |
| Press release before specialist review | Integrity violation |

---

## Project policy

- No automatic generation of “Navier–Stokes solved” claims.
- Do not submit putative solutions directly to CMI before specialist review and peer publication (see `docs/PUBLICATION_STRATEGY.md`).
- Every result must tag **Route A/B/C/D** and evidence level `N*`.
