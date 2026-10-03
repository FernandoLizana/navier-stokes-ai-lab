# Scientific Integrity

This document defines **internal** evidence tags used in this repository.  
**N0–N11 are a project convention**, not an academic standard and not a universal ranking in which every higher index strictly dominates every lower one for all purposes.

## Evidence states (internal tags)

| Level | Meaning (internal) |
|-------|---------------------|
| N0 | Informal idea |
| N1 | Visual observation |
| N2 | Conventional floating-point numerics |
| N3 | Independently reproduced numerics |
| N4 | Arbitrary-precision float |
| N5 | Interval arithmetic |
| N6 | Formal mathematical conjecture (stated) |
| N7 | Hand proof of a lemma |
| N8 | Lean formalization (of the stated finite claim) |
| N9 | Computer-assisted proof with verifiable certificate (as defined here) |
| N10 | Specialist-reviewed |
| N11 | Published and community-accepted |

- Only **N9+** (under this ladder) may be called “computer-assisted proof” **in this project’s language**.
- Only **N11** may be presented as an established mathematical result **in the community sense**.
- No automated claim that continuum Navier–Stokes is solved.

## Accessible evidence labels

Use these plain labels alongside N-tags when writing for newcomers:

| Label | Typical N-range | Means |
|-------|-----------------|--------|
| Idea / conjecture | N0–N6 | Stated claim; not demonstrated |
| Floating-point simulation | N2–N3 | Truncated model run |
| Interval bound (assumptions pending) | N5 | Enclosure under listed hypotheses |
| Certificate structure check | verifier PASS | Software checks only |
| Discrete formal lemma | N8 | Finite / combinatorial statement in Lean |
| Historically reported | — | Written in an old report; re-verify |
| Corrected / invalid / replaced | — | Superseded by repair or audit |
| Independent review documented | N10+ | Only if a review artifact exists |

## Explicit distinctions

- Passing software tests ≠ proving a theorem.
- Verifier `PASS` certifies only the implemented checks and their assumptions.
- A bound too large to close an argument does **not** automatically refute the conjecture.
- A truncated simulation does **not** prove continuum regularity or blow-up.
- Checking finite instances does **not** prove a claim for all IC or all times.
- High precision ≠ rigorous interval enclosure.
- A Lean theorem about declared constants does **not** prove those constants match the analytic object unless that link is proved.
- Presence of Lean files ≠ successful `lake build`. Check `sorry`, axioms, and theorem scope.
- Forced vs unforced equations, domains, and Euler vs Navier–Stokes must stay separated.

## Simulation vs validation vs proof

See also `docs/NUMERICAL_VALIDATION.md` and `docs/COMPUTER_ASSISTED_PROOFS.md`.

## Forbidden confusions (auditor checklist)

Euler vs NS; \(\nu=0\); hyperviscosity; wrong dimension/domain; assuming regularity to prove regularity; ignoring pressure; dropping Fourier tails; truncated blow-up ≠ PDE blow-up; small residual ≠ existence; solution-dependent constants; weak non-uniqueness ≠ smooth blow-up; nonzero force on Route B claims.
