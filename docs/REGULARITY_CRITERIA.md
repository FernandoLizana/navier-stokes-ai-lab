# Regularity Criteria

**Evidence:** N7 for cited classical criteria; NS-MRL does not claim new criteria in Sprint 1.

Continuation / regularity criteria give **conditional** regularity: if a certain quantity remains controlled on \([0,T]\), then the smooth solution extends beyond \(T\). They do **not** by themselves prove Clay Route B.

## 1. Beale–Kato–Majda (1984)

If \(u\) is a smooth solution on \([0,T)\) and

\[
\int_0^T \|\omega(t)\|_{L^\infty}\,dt < \infty,
\]
then the solution extends smoothly past \(T\).

**Implication for numerics:** rapid growth of \(\|\omega\|_\infty\) is a **candidate** signal, not a proof of singularity. Conversely, a solver that “keeps running” does not prove the integral stays finite.

## 2. Prodi–Serrin / Escauriaza–Šeregin–Šverák

If \(u\in L^p(0,T; L^q)\) with \(\frac{2}{p}+\frac{3}{q}=1\), \(q\in(3,\infty]\), then regularity holds. The endpoint \(q=3\) was settled by Escauriaza–Šeregin–Šverák (2003): \(u\in L^\infty(0,T;L^3)\) implies regularity.

## 3. Other classical conditional results (pointers)

- Constantin–Fefferman: geometric conditions on vorticity direction.
- Ladyzhenskaya–Prodi–Serrin family.
- Critical-space smallness ⇒ global regularity (Kato, Koch–Tataru, etc.) — **not** large-data Clay.

## 4. What NS-MRL records as diagnostics (exploratory)

| Quantity | Related criterion | Limitation |
|----------|-------------------|------------|
| \(\|\omega\|_\infty\) | BKM integrand | Grid-dependent; aliasing |
| \(\|\omega\|_{L^2}\) (enstrophy) | Energy cascade proxy | Can grow large without blow-up |
| Spectral slope / analyticity radius proxy | Foias–Temam type heuristics | Not a theorem from float FFT |
| Strain–vorticity alignment | Geometric criteria | Observed ≠ universal |

## 5. Auditor rules

- Do not declare singularity from rapid growth alone.
- Do not declare regularity from successful time-stepping.
- Do not replace \(\nu>0\) by Euler or hyperviscosity when claiming Route B relevance.
- Distinguish Type I vs Type II blow-up **hypotheses** from proofs.

## Primary references

- J. T. Beale, T. Kato, A. Majda, *Comm. Math. Phys.* (1984).
- L. Escauriaza, G. Šeregin, V. Šverák, *Usian. Anal.* (2003).
- J. Serrin, *Arch. Rational Mech. Anal.* (1962).
- P. Constantin, C. Fefferman, *Indiana Univ. Math. J.* (1993).
