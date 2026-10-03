# Theorem Dossier — Finite Galerkin Enstrophy Bounds (L-0001 … L-0006)

**Route:** B (global regularity, T³, exploratory finite-dimensional layer)
**Scope discipline:** Every statement below is about a **finite-dimensional
Galerkin / pseudospectral truncation** on the `N³` grid. None of them is a
statement about the continuum Navier–Stokes PDE, and **none has any Clay
Millennium implication.** They exist to (a) give rigorous ceilings against which
numerical enstrophy-maximization can be audited, and (b) practice the
proof→certificate→verifier pipeline on statements we can actually close.

Evidence ladder tags used here: **N2** = float simulation, **N5** =
interval-certified constant, **N6** = conjecture, **N7** = hand proof for the
truncated system.

---

## Common setup

- Domain `T³ = [0,2π)³`, divergence-free field `u`, energy `E(t)=½‖u‖₂²`,
  enstrophy `Ω(t)=½‖∇u‖₂² = ½‖ω‖₂²`.
- Retained modes: `M` = number of retained wavevectors, `K` = max `|k|`
  (Euclidean) on the retained set. Both are **exact integers/roots** for a
  given grid + dealias + shell.
- Energy is non-increasing (`ν≥0`), so `E(t) ≤ E0` throughout.

---

## L-0001 — Crude cubic bound (N7)

**Statement.** For the finite truncation, dropping viscosity,
`dΩ/dt ≤ C Ω^{3/2}` with `C = 2√2 √(3M)`, hence on `[0,t]`
```
Ω(t) ≤ Ω0 / (1 − ½ C √Ω0 · t)² ,   Ω0 ≤ K²E0,
```
valid while the denominator stays positive. **Weakness:** ignores dissipation;
blows up in finite `t*`.

## L-0002 — Viscous cubic comparison (N7)

Adds the leading viscous sink: `dΩ/dt ≤ −2ν Ω + C Ω^{3/2}`. Comparison ODE can
still diverge in finite time when `a√Ω0 > ν`; **auditor note:** that divergence
is an *estimate artifact of the truncated ODE bound*, **not** a physical blow-up.

## L-0003 — Uniform-in-time bound (N7)

Uses the enstrophy–energy Cauchy inequality `‖∇ω‖₂² ≥ 2Ω²/E`:
```
dΩ/dt ≤ −ν (2Ω²/E0) + C Ω^{3/2}   ⟹   Ω(t) ≤ max( Ω0 , 6 M E0² / ν² ).
```
First **globally bounded** ceiling of the family (attracting equilibrium). Not
uniform as `N→∞` (`M→∞`).

## L-0004 — Spectral-support refinement (N7 shell-restricted / N2 diagnostic)

Replaces `(K,M)` by shell-effective `(K_eff, M_eff)` when the field's Fourier
support is concentrated. Shell-restricted version is N7; the field-conditional
version (measured support of an actual trajectory) is only an **N2 diagnostic**.

## L-0005 — Hard Galerkin stays in shell (N7)

If the dynamics are **truncated to** `|k| ≤ K0` at every step (projection inside
the timestep), the support never leaves the shell. This is the hypothesis that
makes shell bounds airtight — and is precisely **not** what C-S-000x assume
(their ICs are shelled but evolution runs on the full dealiased grid).

## L-0006 — Improved short-time shell bound (N7 + N5 constant)

For hard shell truncation, `‖∇u‖_∞ ≤ √(3M)√(2Ω)` gives `dΩ/dt ≤ C Ω^{3/2}` and
the algebraic envelope
```
Ω(t) ≤ Ω0 / (1 − a√Ω0 · t)² ,   a = √2 √(3M),   Ω0 = K²E0,
```
while `a√Ω0 t < 1`; otherwise falls back to L-0001/L-0003. **Certified instance:**
see below.

---

## Certificate CERT-L0006-shell-k2 (N5)

Machine-checked **interval-arithmetic** enclosure of the L-0006 algebraic bound
for the hard shell `|k| ≤ 2` on `N=12`, `E0=0.5`, `t∈[0,0.02]`:

| quantity | value |
|---|---|
| `K²` (exact) | 4 |
| `M` (exact) | 33 |
| `a√Ω0·t` upper bound | `0.39799…` (< 1 ⇒ closes) |
| **guaranteed** `Ω(t) ≤` | **5.5186…** |

- Generator: `python/ns_exploration/validation/l0006_certificate.py`
- Independent verifier: `python/ns_exploration/validation/l0006_certificate_verify.py`
- Stored artifact: `certificates/CERT-L0006-shell-k2.json`.

## Certificate CERT-L0003-shell-k3 (N5) — fallback

For `|k| ≤ 3` (`K²=9`, `M=123`), the L-0006 algebraic product
`a√Ω0·t ≈ 1.153 > 1` **fails to close** at `t=0.02`. The uniform L-0003
ceiling is used instead:

| quantity | value |
|---|---|
| `K²`, `M` (exact) | 9, 123 |
| `Ω0 ≤ K²E0` (hi) | 4.5 |
| `Ω_eq = 6ME0²/ν²` (hi) | 18450 |
| **guaranteed** `Ω(t) ≤` (all `t≥0`) | **18450** |

- Generator / verifier: `validation/l0003_certificate*.py`
- Artifact: `certificates/CERT-L0003-shell-k3.json`

## Certificate CERT-L0003-full-dealias-N12 (N5)

Same L-0003 formula on the **full 2/3 dealias mask** (standard pseudospectral
Galerkin support), not a Euclidean shell:

| quantity | value |
|---|---|
| `K²`, `M` (exact) | 27, 343 |
| `Ω0 ≤ K²E0` (hi) | 13.5 |
| `Ω_eq` (hi) | 51450 |
| **guaranteed** `Ω(t) ≤` (all `t≥0`) | **51450** |

- Artifact: `certificates/CERT-L0003-full-dealias-N12.json`

## Certificate CERT-L0003-full-dealias-N16 (N5) — growth witness

| quantity | value |
|---|---|
| `K²`, `M` (exact) | 75, 1331 |
| `Ω0_hi` | 37.5 |
| **guaranteed** `Ω(t) ≤` | **199650** |
| growth vs N12 | `M₁₆/M₁₂ = 1331/343 ≈ 3.88 = cap₁₆/cap₁₂` |

Documents the expected `∼ M/ν²` blow-up of the uniform Galerkin constant as
resolution increases. Still finite-dimensional only; **no Clay claim**.
Artifact: `certificates/CERT-L0003-full-dealias-N16.json`.

## Lean skeleton (toward N8)

Package `lean/NSGalerkin` (`lake build`):

| module | statement |
|---|---|
| `EnstrophyEnergy` | Nat: `Σ k²·e ≤ K² · Σ e` |
| `EnstrophyEnergyRat` | ℚ≥0 with common den `D`: same inequality for `eᵢ = nᵢ/D` (Fourier `|û|²≥0`) |
| `ComparisonODE` | sign of comparison force; `max(Ω0,Ω_eq)` cap; product form |
| `Certificates` | exact ℚ constants: shell-k3=18450, N12=51450, **N16=199650**; growth `N12 < N16` |

This is an **N8 skeleton for elementary identities**, not a formalization of
NS or Clay. The differential inequality itself remains hand N7. The
`Certificates` module cross-verifies, inside Lean, the same numeric constants
that the Python interval verifier guards; `test_cert_lean_bridge.py` pins the
Python certificates to those Lean-proved integers so drift on either side
fails a test.

---

## Relationship to the conjecture line

- **C-0001** (all-IC enstrophy bound) — refuted (N2), archived.
- **C-0002** (all-IC, M≈20.32) — **refuted** (N2) by Stokes `|k|²=75` mode.
- **C-0003** (all-IC, M=37.5) — **proved** by L-0018.
- **C-0004** (all-IC, N≤24, M=73.5) — **proved** by L-0018.
- **C-0005** (all-IC, N≤24, M≈61.24 at T=0.02) — **proved** by L-0024.
- **C-0006** (all-IC, N≤32, M≈67.77 at T=0.02) — **proved** by L-0025.
- **C-0007** (all-IC, N≤24, M≈41.28) — **active**; subclasses C-R-0002…0012; techos L-0026…L-0048 (L-0048: fullsym Shor).
- **C-R-0001** (mono-radial IC, N≤24, M≈71.73) — **proved** by L-0020+L-0021.
- **C-R-0002** (Ω0≥Ω★≈18.87, N≤24, M≈41.28) — **proved** by L-0024+L-0026.
- **C-R-0003** (mono-radial α_r=0 shells, M≈41.28) — **proved** (Stokes–Duhamel N_*=0).
- **C-R-0004** (Fourier support on shells {1,2,3}, M≈41.28) — **proved** by L-0038+L-0026+L-0027.
- **C-R-0005** (14-shell support A incl. {1,2,3,147,…}, M≈41.28) — **proved** by L-0039+L-0026+L-0027.
- **C-R-0006** (Fourier support on shells {1..5}, M≈41.28) — **proved** by L-0040 Sym Shor + L-0026+L-0027.
- **C-R-0007** (8-shell support incl. r=6, M≈41.28) — **proved** by L-0041 Sym Shor + L-0026+L-0027.
- **C-R-0008** (one-pol shells {1,2,3,4,5,6,8}, M≈41.28) — **proved** by L-0043+L-0026+L-0027.
- **C-R-0009** (Hermitian shells {1..6}, M≈41.28) — **proved** by L-0044 physical fullsym + L-0026+L-0027.
- **C-R-0010** (Hermitian shells \|k\|²≤25, M≈41.28) — **proved** by L-0045 streaming fullsym + L-0026+L-0027.
- **C-R-0011** (30 Hermitian shells incl. 147, M≈41.28) — **proved** by L-0046 greedy fullsym + L-0026+L-0027.
- **C-R-0012** (41 Hermitian shells incl. 147, M≈41.28) — **proved** by L-0047 greedy fullsym + L-0026+L-0027.
- **C-S-0001** (shell IC, full-grid evolution) — refuted (N2), archived.
- **C-S-0002** (shell IC, `M≈48.16`) — **proved** by L-0018 (envelope 37.5).

The lemmas L-0001…L-0006 bound the **hard-truncated** system; the C-S conjectures
concern the **full-grid** evolution of shelled ICs. **L-0007** is the bridge:
same full-grid stretch constants, but initial ceiling `Ω0 ≤ K_IC² E0`.
**L-0018** supersedes the need for growth ODEs when `M ≥ K²(N) E0`.

## L-0007 — Shell IC → full dealias grid (N7 + N5)

For ICs with `|k|_∞ ≤ k_inf` (or Euclidean `|k| ≤ k_eucl`) evolving on the
full 2/3 dealias grid:

```
dΩ/dt ≤ 2 a Ω^{3/2},   a = √2 √(3 M_full),   Ω(0) ≤ K_IC² E0,
```

plus exponential fallback `Ω(t) ≤ Ω0 e^{α t}` with
`α = 2 K_full √(3 M_full) √(2 E0)`, and the L-0003 full-grid uniform ceiling.

| case | best proved cap | vs L-0003-full | vs C-S-0002 M≈48 |
|---|---|---|---|
| N=12, eucl. k≤2 | **≈1573** (exp) | 32× sharper | still ≫ 48 |
| N=12, ℓ^∞ k≤4 | ≈1.89e4 (exp) | sharper | ≫ 48 |
| N=16, ℓ^∞ k≤4 (C-S class) | **199650** (L-0003) | same | gap ≈ **4146×** |

N5 artifact: `certificates/CERT-L0007-exp-euclidean-k2-N12.json` (bound ≤ 1572.53).

**Honesty:** L-0007 does **not** prove C-S-0002; it makes the gap quantitative.

## L-0008 — Cascade via `E_H(0)=0` (N7 + N5)

High modes start empty. With `B = K_IC² E0` and
`A = √2 √(3 M_full) √E0`,
```
Ω(t) ≤ B · cosh²(A K_full t).
```

| case | L-0008 cosh | L-0007 | improvement |
|---|---|---|---|
| N=12, eucl. k≤2 | **≈394** | ≈1573 | ~4× |
| N=12, ℓ^∞ k≤4 | ≈4730 | ≈1.89e4 | ~4× |
| N=16, ℓ^∞ k≤4 | explodes | **199650** (L-0003) | cosh loses; gap to M≈48 still ~4146× |

N5: `certificates/CERT-L0008-cosh-euclidean-k2-N12.json`.

Still does **not** prove C-S-0002. The remaining gap is the crude `√(3M_full)`
embedding on the whole grid once any high mode exists.

## L-0009 — Two-scale embedding (N7 + N5)

```
||u||_∞ ≤ √(3 M_L)√(2E0) + √(3 M_H)√(2 E_H) = U_L + W √E_H,
```
then `z' ≤ (U_L+Wz)√(B+K²z²)` with closed form via `z=(√B/K)sinh u`:
`Ω = B cosh²(u(t))`.

| case | L-0009 | L-0008 | notes |
|---|---|---|---|
| N=12, eucl. k≤2 | **≈103.58** | ≈394 | ~3.8× sharper; N5 certified |
| N=16, ℓ∞≤4 (C-S) | does not close | 199650 | ODE saturates before T=0.02 |

Artifact: `certificates/CERT-L0009-twoscale-euclidean-k2-N12.json`.
C-S-0002 remains unproved (gap still ~4146× at N=16).

## L-0010 — Duhamel short-time bootstrap (N7 + N5)

Bootstrap `Ω ≤ R` on `[0,T]` via
`E_H ≤ T² (U_L + W√ε)² R` with `ε=(R-B)/K_full²`.

| statement | status |
|---|---|
| C-S-0002 at T=0.02, R≈48.16 | **still open** |
| **C-S-0002-SHORT**: same R on `[0, T*]`, `T*≈0.001011` | **proved restricted (N7)** |
| N5 cert at `T_cert=0.99 T*` | `CERT-L0010-duhamel-CS0002-short-N16` |

This is the first *positive* theorem in the C-S line: the conjectured constant
holds on a shorter explicit horizon. Extending `T*` → 0.02 remains open.

## L-0011 — Sharpened short-time (div-free + viscous)

| quantity | L-0010 | L-0011 |
|---|---|---|
| rigorous `T*` for R≈48.16 | 0.001011 | **0.001240** (~1.23×) |
| conditional L×L `T*` (N6) | — | 0.003046 |
| conditional `Ω(0.02)` L×L (N6) | — | ≈1023 (≪ 199650) |

N5: `CERT-L0011-duhamel-CS0002-short-N16`. C-S-0002 at T=0.02 still open.

## L-0012 — NS cancellations for high-mode energy (N7 + N5)

After ⟨u_H,(u_L·∇)u_H⟩=⟨u_H,(u_H·∇)u_H⟩=0,
```
dE_H/dt ≤ -2νκ_H E_H + 2 U_L√B √E_H + 2 W√B E_H
z' ≤ (W√B − νκ_H) z + U_L√B,   z=√E_H
Ω ≤ B + K_full² z²
```

| quantity | L-0011 | L-0012 |
|---|---|---|
| rigorous `T*` for R≈48.16 | 0.001240 | **0.002284** (~1.84×) |
| fraction of parent T=0.02 | ~6.2% | ~11.4% |

N5: `CERT-L0012-cancel-CS0002-short-N16`. Parent C-S-0002 at T=0.02 still open
(structural: full high-mode ‖u‖_∞ contribution cannot close to 0.02 for finite R).

## L-0013 — Vector Fourier CS embedding (N7 + N5)

Corrects the L-0011/12 factors `√(2M)√(2E)` / `√(4 M_H)` to the sharp
Cauchy–Schwarz mode sum
`‖u‖_∞ ≤ √M √(2E)` (exclude `k=0`), so `U_L=√M_L√(2E0)`, `W=√(2 M_H)`.

| quantity | L-0012 | L-0013 |
|---|---|---|
| `U_L` | ≈38.18 | **≈26.98** |
| `W` | ≈49.07 | **≈34.70** |
| rigorous `T*` for R≈48.16 | 0.002284 | **≈0.00324** (~1.42×) |
| fraction of parent T=0.02 | ~11.4% | ~16.2% |

N5: `CERT-L0013-embed-CS0002-short-N16`. Parent C-S-0002 at T=0.02 still open
(even L×L-only majorant needs `U_L ≲ 5.94` to close).

## L-0014 — Shell→high triad multiplicity (N7 + N5)

Convolution counting: `R_H = max_{k∈H} #{(p,q)∈S_L² : p+q=k}` (=324 for N=16, ℓ∞≤4).
Then `‖P_H(u_L·∇)u_L‖_ℓ₂ ≤ 2√(R_H E Ω)`, so the L-0012 ODE runs with
`U_eff = √R_H √(2E0)` (=18) instead of `√M_L √(2E0)` (=26.98).

| quantity | L-0013 | L-0014 |
|---|---|---|
| production `U` | ≈26.98 | **18** (`R_H=324`) |
| rigorous `T*` for R≈48.16 | 0.003235 | **≈0.00437** (~1.35×) |
| fraction of parent T=0.02 | ~16.2% | ~21.8% |

N5: `CERT-L0014-triad-CS0002-short-N16`. Parent C-S-0002 at T=0.02 still open
(need `U_eff ≲ 5.94`; empirical bilinear ratios suggest further room).

## L-0015 — S_max energy triad bound (N7 + N5)

Row-weighted energy estimate:
`‖P_H N‖_ℓ₂ ≤ 2 E √S_max` with
`S_max = max_p Σ_{q: p+q∈H} |q|²` (=6132 for N=16, ℓ∞≤4).
Gives `U_eff≈11.30` (vs L-0014's 18 from `R_H`).

| quantity | L-0014 | L-0015 |
|---|---|---|
| production `U_eff` | 18 | **≈11.30** |
| rigorous `T*` for R≈48.16 | 0.004367 | **≈0.00597** (~1.37×) |
| fraction of parent T=0.02 | ~21.8% | ~29.8% |

N5: `CERT-L0015-smax-CS0002-short-N16`. Parent C-S-0002 at T=0.02 still open
(empirical `‖N‖/√(EΩ)≲0.7` vs certified scale ~`N_max/√(E0 B)≈22.6`).

## L-0016 — Radial Young bound (N7 + N5)

`‖A∗B‖_ℓ₂ ≤ ‖A‖₂‖B‖₁` with radial grouping:
`ρ_★ = max_r (r m_r) = 1392` (r=29, m_r=48) ⇒ `N_max = 2 E0 √ρ_★ ≈ 37.31`,
`U_eff ≈ 5.39` (vs L-0015 ≈11.30).

| quantity | L-0015 | L-0016 |
|---|---|---|
| production `U_eff` | ≈11.30 | **≈5.39** |
| rigorous `T*` for R≈48.16 | 0.005967 | **≈0.00912** (~1.53×) |
| fraction of parent T=0.02 | ~29.8% | ~45.6% |

N5: `CERT-L0016-radial-CS0002-short-N16`. Parent C-S-0002 at T=0.02 still open
(with current `W`, need `U_eff≲0.5` to cross 0.02; L×L-only needs `≲5.94`).

## L-0017 — Radial Young for H×L cross (N7 + N5)

Replaces `W√B` in the cancellation ODE by
`γ_cross = √(2 ρ_★ E0) ≈ 37.31` (≪ `W√B ≈ 170`), via Young on `(u_H·∇)u_L`.

| quantity | L-0016 | L-0017 |
|---|---|---|
| cross coeff in `z'` | ≈170 (`W√B`) | **≈37.31** (`γ_cross`) |
| rigorous `T*` for R≈48.16 | 0.009116 | **≈0.01606** (~1.76×) |
| fraction of parent T=0.02 | ~45.6% | **~80.3%** |
| `Ω(0.02)` comparison | ≈1431 | ≈67.6 |

N5: `CERT-L0017-cross-CS0002-short-N16`. Parent C-S-0002 at T=0.02 still open
*under the L-0012–17 comparison ODE* (loose `Ω ≤ K_IC² E0 + K_full² E_H`).

## L-0018 — Spectral enstrophy envelope (N7 + N5) → **C-S-0002 proved**

Parseval on the dealias mask: `Ω ≤ K²(N) E`. With `E(t)≤E0` and
`K²(16)=75` (worst among `N≤16`),
`Ω(t) ≤ 37.5 < M≈48.16` for all `t≥0`. Shell IC is a subclass.

| quantity | L-0017 comparison @0.02 | **L-0018 envelope** |
|---|---|---|
| `Ω` upper bound | ≈67.6 (loose ODE) | **≤37.5** |
| closes C-S-0002 | no | **yes (all t≥0)** |

N5: `CERT-L0018-envelope-CS0002-N16`. Does **not** close the *old* `C-0002` (`M≈20.32`);
that conjecture is **refuted** by a Stokes `|k|²=75` mode (`Ω(0.02)≈27.78`).
Successor **C-0003** takes `M=37.5` and is proved by the same envelope
(`CERT-L0018-envelope-C0003-N16`). Not continuum / not Clay.

## L-0018 @ N≤24 → **C-0004 proved**

Same envelope with `K²(24)=147`: `Ω(t) ≤ 73.5` for all `t≥0` on N≤24.
N5: `CERT-L0018-envelope-C0004-N24`.

## L-0019 — Stokes spectral majorant (linear)

`S(t)=max_k |k|² e^{-2ν|k|² t}`; `Ω_Stokes(t)≤ S(t) E0`.
At N=24, T=0.02: floor ≈40.825 (attained). Places **C-0005** (`M≈61.24`)
in the open gap `(floor, envelope)=(40.825, 73.5)`. Does **not** bound nonlinear NS.
N5: `CERT-L0019-stokes-floor-N24`.

## L-0020 — Stokes–Duhamel H¹ majorant

\[
\Omega(T)\le \tfrac12\bigl(\sqrt{2 S(T)E_0}+N_* I_\sigma\bigr)^2.
\]

Young `N_*=2 E0√ρ_★≈78.7` ⇒ `Ω_H1≈327` (worse than envelope).  
**C-0005 closes if `N_*≤N_crit≈9.666`** (empirical `‖N‖≲2`).  
N5: `CERT-L0020-duhamel-H1-N24`.

## L-0021 — Mono-radial quartic Gram

On each radial shell, Shor-style Frobenius relaxation yields `α_r ≥ max ‖N‖`
at `‖u‖_ℓ₂=1`. At N=24: `α_★ = max_r α_r = 14` (r=74). With L-0020:
`Ω_H1≈71.73 < 73.5` (beats envelope) but `>61.24` (does not close C-0005).

**C-R-0001** (mono-radial ICs only): proved at `M≈71.73`.  
N5: `CERT-L0021-quartic-radial-N24`. (C-0005 later closed by L-0024.)

## L-0022 — Gamma / N_* bootstrap arithmetic

Conditional: `‖N‖ ≤ γ √(E Ω)` and bootstrap `Ω≤M` ⇒ L-0020 closes C-0005 iff
`γ ≤ γ_crit ≈ 1.7469`. Adversarial N2: `‖N‖ ≳ 35 ≫ N_crit`, `γ ≳ 10 ≫ γ_crit`,
while `Ω(T) ≲ 41 ≪ M`. **Uniform** N_*/γ cannot close all-IC C-0005.
N5: `CERT-L0022-gamma-bootstrap-N24`.

## L-0023 — Cubic-dissipation ODE (path to C-0005)

From `Σ|k|⁴ ‖û‖² ≥ 2Ω²/E` and `E≤E0`,
`dΩ/dt ≤ -(2ν/E0) Ω² + ⟨ω, curl N⟩`.
If `|⟨ω, curl N⟩| ≤ C Ω^{3/2}` and `C ≤ C_★ ≈ 2.154`, then
`Ω(0.02) ≤ M_{C-0005}` from `Ω0 ≤ 73.5`. Empirics `C_emp ≪ C_★`.
Proved `C≤C_★` open (L-0002's `C~√M` too weak). N5: `CERT-L0023-cubic-dissipation-N24`.
Superseded for C-0005 closure by **L-0024**.

## L-0024 — Spectral defect ⇒ C-0005 proved

Max shell `|k|²=147` (8 corner modes) has no dealias self-triads ⇒ `N(u_K)=0`.
Defect `δ=K²E-Ω ≥ 13 E_<` and
`‖N‖ ≤ a(E)√δ + b δ` feed a coupled `(E,Ω)` comparison ODE.
Worst-case `Ω(0.02) ≤ ≈41.74 < M≈61.24`. N5: `CERT-L0024-spectral-defect-C0005-N24`.

## L-0025 — Multi-N spectral defect ⇒ C-0006

Same defect ODE at N=32 (`K²=300`, gap=19, no self-triads): worst `Ω(0.02)≈46.15`
closes **C-0006** at `M=1.5×Stokes_floor≈67.77`. N5: `CERT-L0025-spectral-defect-C0006-N32`.
Remaining open gap on N=24: **C-0007** (`M≈41.28` between floor and L-0024 majorant).
**L-0026** shows the high-Ω0 slab (`Ω0≥Ω★≈18.87`) already meets M via L-0024; the low-Ω0
slab is a structural techo for defect-ODE sharpenings (C-0007 still open).
High slab packaged as **C-R-0002**.

## L-0027 — Low-slab cubic conditional

On `Ω0≤Ω★`, the L-0023 ODE closes `Ω(T)≤M` if stretch obeys `|⟨ω,curl N⟩|≤C Ω^{3/2}`
with `C≤C_†≈9.562` (emp ≪ C_†; proved C open). All-IC cubic at `C=0` already exceeds M,
so the high slab must stay on L-0024. N5: `CERT-L0027-low-slab-cubic-C0007-N24`.

## L-0028 — Embedding techo for C_†

L-0002-style `C=2√2√(3M)` meets `C_†` only for `M≤M_★≈3.81`. The `|k|²=1` shell already
has `M=6` (`C≈12.0`); `|k|²≤2` gives `C≈20.8`. Embeddings cannot close the low slab —
need a triad-tensor stretch bound. N5: `CERT-L0028-embedding-techo-C0007-N24`.

## L-0029 — Frobenius Jacobian embedding

`|stretch|≤2Ω‖∇u‖_{∞,F}` with `‖∇u‖_{∞,F}≤√M√(2Ω)` gives `C_F=2√(2M)` (factor `√3`
below L-0002). Full mask: `C_F≈164.3≫C_†`. Needs `M≤11.43`; r=1 meets (`C_F≈6.93`) but
`λ_max=1≪λ★≈37.75`. Still no all-IC low-slab closure. N5: `CERT-L0029-frobenius-stretch-C0007-N24`.

## L-0030 — Full-mask triad ‖N‖

`R_★=max_k #{p+q=k}=3148` on nonzero dealias modes ⇒ `‖N‖≤2√(R_★ E Ω)`
(beats Young at Ω∼E: ≈56.1 vs ≈78.7). Absolute cubic majorant gives `C≲158.7`, still
`≫C_†`; low-slab hybrid ODE floors ≳144. Need weighted/SOS triad bounds.
N5: `CERT-L0030-triad-N-C0007-N24`.

## L-0031 / L-0032 / L-0033 — Exhaustive techos

- **L-0031**: compatible-S defect majorant → worst `Ω(T)≈41.825` (worse than L-0024).
- **L-0032**: full-mask `S_max≈1.70×10⁵` ⇒ `N≲412` ≫ Young.
- **L-0033**: absolute weighted cubic Frobenius ⇒ `C≲8992` ≫ `C_†`.
- **C-R-0003**: mono-radial shells with `α_r=0` (incl. 147) meet M via Stokes–Duhamel.

N5: `CERT-L0031-0033-techos-C0007-N24`. Remaining: signed/SOS stretch tensor.

## L-0034 — Shell-factorized ‖N‖

`‖N‖≤2 A_★(E,Ω) B_★(E,Ω)` with two-point maxima of `Σ√(m e)` and `Σ√(r e)`.
At `Ω=E=0.5`: `N≲2.45≪Young≈78.7`, but `C_eff≳84≫C_†` and low-slab ODE `≳87>M`.
N5: `CERT-L0034-shell-N-C0007-N24`.

## L-0035 — Shell two-point ‖∇u‖_{∞,F}

`‖∇u‖_{∞,F}≤Σ√(r m_r)√(2 e_r)` maximized at two-point spectra ⇒ `C=2 g_∞/√Ω`.
At `Ω=E`: `C≈6.93≤C_†`, but on the low slab `C_max≳36` and `Ω_c≈E` only.
Pure ginf ODE ≳80; hybrid min(ginf, defect, shell-N, Young) still ≳55 > M.
N5: `CERT-L0035-ginf-stretch-C0007-N24`. Remaining: signed/SOS triad stretch.

## L-0036 — Signed Shor / matricization stretch

Ω-weighted polarization cubic `T(z,z,z)=z·L(zzᵀ)` ⇒ `C≤C_Shor=2√2‖L‖`.
Power-iteration lower bounds: `C_Shor≳42.8` (N=12), `≳141.8` (N=24), both `≫C_†`.
Shor relaxation is too loose; need SOS / rank-constrained or shell-diff blocks.
N5: `CERT-L0036-signed-shor-C0007-N12`.

## L-0037 — Shell-difference / bi-shell Shor

Mono-radial: `stretch = r Σ T_k = 0` exactly. Bi-shell Shor on N=24: near pairs
can meet `C_†` (e.g. `(1,2)≲6.7`), but sample max `≳32` (e.g. `(2,134)`) `≫C_†`.
Shell-diff is real but not uniform enough for all-IC. Next: multi-shell SOS/rank.
N5: `CERT-L0037-shell-diff-C0007-N24`.

## L-0038 — Exact SVD Shor on low-shell bands (+ C-R-0004)

Dense SVD of the Ω-weighted stretch matricization on shell bands:
`{1,2}→6.633`, `{1,2,3}→9.165`, `{1,2,5}→9.550` all `≤C_†`; `{1,2,3,4}→10.066>C_†`.
**C-R-0004** (support on `{1,2,3}`): C-0007 proved via L-0038+L-0026+L-0027.
All-IC still open. N5: `CERT-L0038-band-svd-C0007-N24`.

## L-0039 — Sparse SVD multi-shell supports (+ C-R-0005)

Sparse SVD of stretch matricization on greedy enlargements of L-0038 bands:
support A (14 shells, incl. `{1,2,3,147}`) `C_Shor≈9.519≤C_†`; support B
(14 shells from `{1,3,4,6}`) `≈9.548≤C_†`; quad `{1,3,4,6}≈8.485`.
**C-R-0005** on support A. All-IC still open. N5: `CERT-L0039-sparse-support-C0007-N24`.

## L-0040 — Sym-restricted Shor (+ C-R-0006)

Restrict matricization to symmetric `Z=zzᵀ`: `C_Sym=2√2‖L|_{Sym}‖`.
`{1..5}→≈9.015≤C_†`; `{1..6}→≈11.34>C_†`. Unrestricted Shor failed already on `{1..4}`.
**C-R-0006**: support on shells `{1,2,3,4,5}`. N5: `CERT-L0040-sym-shor-C0007-N24`.

## L-0041 — Shell-6 Sym catalog (+ C-R-0007)

Consecutive `{1..6}` has `C_Sym≈11.34>C_†` (techo; rank-1 still ≪C_†). Leave-one-out
seeds `{1,2,3,4,6}`, `{2,3,4,5,6}`, `{1,3,4,5,6}` meet C_†; greedy from the first
gives **C-R-0007** on `{1,2,3,4,6,11,12,25}` (`C_Sym≈9.414`). N5: `CERT-L0041-shell6-sym-C0007-N24`.

## L-0042 — Sym-strengthening techo on {1..6}

Tried to beat `C_Sym({1..6})≈11.34`: column-triangle `≈80`, A/B hybrid `≈16.2`
(both worse); rank-1 lower `≈0.23≪C_†`. Sym remains best proved majorant and still
misses C_†. Need SOS / true rank-constrained upper bound. N5: `CERT-L0042-sym-strengthen-techo-C0007-N24`.

## L-0043 — One-pol Sym (+ C-R-0008) + shell-block/Ky-Fan techos

Single frozen polarization per wavevector (first `_pol_basis` branch) yields
`C_Sym≈8.476≤C_†` on shells `{1,2,3,4,5,6,8}` (**C-R-0008**); adding shell 9
exceeds C_† (`≈9.849`). Two-pol shell-block majorant lo `≳14.5` and 4-term Ky-Fan
`≳12.6` both fail to close two-pol `{1..6}`. N5: `CERT-L0043-onepol-shellblock-C0007-N24`.

## L-0044 — Physical Hermitian fullsym (+ C-R-0009)

Audit (Sprint A/B): the physical stretch uses the triad factor `i` and a real
Hermitian `(c,s)` basis (`‖z‖²=2Ω`). FFT `stretch_inner` matches the cubic tensor
to `≲1e-13`. Full symmetrization then gives
`C_fullsym({1..6})≈2.993≤C_†` while `(i,j)`-only Sym remains `≈11.344>C_†`.
**C-R-0009** on Hermitian support `{1,2,3,4,5,6}`. Exploratory: `C_fullsym` still
`≤C_†` through shells `r≤20` (`≈7.70`); dense `G` OOMs beyond. All-IC open.
N5: `CERT-L0044-physical-fullsym-C0007-N24`.

## L-0045 — Streaming physical fullsym (+ C-R-0010)

Stream `M=flatten(full_sym(G))` without storing `G`. Maximal consecutive band:
`r_max=25` → `C_fullsym≈9.119≤C_†` (**C-R-0010**, 22 shells, D=1028);
`r_max=26` → `≈9.779>C_†` (techo). Full dealias (`D≈6748`) still open.
N5: `CERT-L0045-streaming-fullsym-C0007-N24`.

## L-0046 — Greedy high-shell fullsym (+ C-R-0011)

Greedy add smallest high shells to C-R-0010 under an 8GB `M` cap. Witness
30 shells incl. `{32,37,40,48,72,98,108,147}` with `C_fullsym≈9.317≤C_†`.
**C-R-0011**. All-IC open. N5: `CERT-L0046-greedy-fullsym-C0007-N24`.

## L-0047 — Extended greedy fullsym (+ C-R-0012)

Batch-greedy under 12GB float32 `M` + Gram `σmax`. Witness 41 shells
(C-R-0011 ∪ `{43,44,52,58,67,68,73,76,88,102,114}`) with `C_fullsym≈9.504≤C_†`.
Rejected `{85,97,107}`; further shells skipped by RAM at `D≈1772`.
**C-R-0012**. All-IC open. N5: `CERT-L0047-greedy-fullsym-C0007-N24`.

## L-0048 — Full-dealias sparse fullsym Shor techo

Sparse CSR assembly of physical fullsym `M` on all dealias shells
(`D=6748`, `nnz≈3.78×10⁷`) + Gram `σmax` gives `C_fullsym≈25.926≫C_†`.
Physical fullsym Shor **cannot** certify all-IC low slab (method techo;
does not refute C-0007). Validated vs dense on `{1..6}` (rel≈1e-16).

## L-0049 — Rank-1 lower bound on full-dealias cubic

Tensor power on sparse `M` (64 starts): `C_rank1≈3.446` (FFT check ≈3.23).
`C_rank1 ≪ C_† ≪ C_fullsym` with gap `C_fullsym/C_rank1≈7.5`. Does **not**
refute C-0007; shows Shor is very loose — next attack needs a tighter
certificate (SOS / hierarchical / structure-exploiting).
