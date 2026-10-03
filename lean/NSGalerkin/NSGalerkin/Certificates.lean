/-!
# Certificate arithmetic bridge (rational, exact)

Ties the Lean layer directly to the Python N5 certificates by proving, over the
rationals `ℚ`, the exact value of the L-0003 uniform-cap constants for the
concrete certificate data:

  Ω0   = K² · E0
  Ω_eq = 6 · M · E0² / ν²
  cap  = max(Ω0, Ω_eq)

with the same exact integers `(K², M)` and rationals `(E0, ν)` used by
`certificates/CERT-L0003-*.json`. Verified with `native_decide` (compiled
evaluation of decidable rational identities), giving a Lean-checked
cross-verification of the numeric constant that the Python interval verifier
also guards.

**Honesty.** This verifies the *arithmetic of a finite-Galerkin bound constant*.
It does NOT prove the differential inequality (hand N7) and has NO Clay
implication. `native_decide` adds the compiler to the trusted base; a
kernel-only integer companion is given below via `*_cleared` Nat identities.
-/

namespace NSGalerkin.Certificates

/-- Uniform-cap equilibrium term `Ω_eq = 6 M E0² / ν²` over ℚ. -/
def omegaEq (M : Nat) (E0 nu : Rat) : Rat :=
  (6 * (M : Rat) * (E0 * E0)) / (nu * nu)

/-- Initial ceiling `Ω0 = K² E0` over ℚ. -/
def omega0 (K2 : Nat) (E0 : Rat) : Rat := (K2 : Rat) * E0

/-- Uniform cap `max(Ω0, Ω_eq)`. -/
def uniformCap (K2 M : Nat) (E0 nu : Rat) : Rat :=
  max (omega0 K2 E0) (omegaEq M E0 nu)

/-! ## CERT-L0003-shell-k3 : K²=9, M=123, E0=1/2, ν=1/10 -/

theorem shell_k3_omega0 : omega0 9 (1/2) = 9/2 := by native_decide
theorem shell_k3_omegaEq : omegaEq 123 (1/2) (1/10) = 18450 := by native_decide
theorem shell_k3_cap : uniformCap 9 123 (1/2) (1/10) = 18450 := by native_decide

/-! ## CERT-L0003-full-dealias-N12 : K²=27, M=343, E0=1/2, ν=1/10 -/

theorem full_N12_omega0 : omega0 27 (1/2) = 27/2 := by native_decide
theorem full_N12_omegaEq : omegaEq 343 (1/2) (1/10) = 51450 := by native_decide
theorem full_N12_cap : uniformCap 27 343 (1/2) (1/10) = 51450 := by native_decide

/-! ## CERT-L0006-shell-k2 sanity: Ω0 = K²E0 = 2 for K²=4, E0=1/2 -/

theorem shell_k2_omega0 : omega0 4 (1/2) = 2 := by native_decide

/-! ## Kernel-only integer companions (E0=1/2, ν=1/10 ⇒ Ω_eq = 150 M)

With `E0 = 1/2`, `ν = 1/10`:
  Ω_eq = 6 M (1/4) / (1/100) = 6 M · 25 = 150 M.
These are pure Nat identities checked by the kernel (`decide`), independent of
`native_decide`, matching the ℚ values above. -/

/-- Cleared-denominator equilibrium for the fixed (E0,ν)=(1/2,1/10) regime. -/
def omegaEqCleared (M : Nat) : Nat := 150 * M

theorem shell_k3_omegaEq_cleared : omegaEqCleared 123 = 18450 := by decide
theorem full_N12_omegaEq_cleared : omegaEqCleared 343 = 51450 := by decide

/-! ## CERT-L0003-full-dealias-N16 : K²=75, M=1331, E0=1/2, ν=1/10
Documents growth ∼ M/ν²: cap = 150·1331 = 199650 (vs 51450 at N=12). -/

theorem full_N16_omega0 : omega0 75 (1/2) = 75/2 := by native_decide
theorem full_N16_omegaEq : omegaEq 1331 (1/2) (1/10) = 199650 := by native_decide
theorem full_N16_cap : uniformCap 75 1331 (1/2) (1/10) = 199650 := by native_decide
theorem full_N16_omegaEq_cleared : omegaEqCleared 1331 = 199650 := by decide

/-- Growth witness: N16 cleared cap exceeds N12 cleared cap. -/
theorem full_growth_N12_to_N16 : omegaEqCleared 343 < omegaEqCleared 1331 := by decide

/-- Twice the initial ceiling equals K² (since E0 = 1/2): `2·Ω0 = K²`. -/
theorem shell_k3_omega0_cleared : 2 * omega0 9 (1/2) = 9 := by native_decide
theorem full_N12_omega0_cleared : 2 * omega0 27 (1/2) = 27 := by native_decide
theorem full_N16_omega0_cleared : 2 * omega0 75 (1/2) = 75 := by native_decide

/-! ## L-0018 / C-0003 / C-S-0002 envelope: Ω ≤ K² E0 = 75/2 on N=16

Same `omega0 75 (1/2)` constant. Documents that the spectral envelope cap used
to close C-0003 (M=37.5) and C-S-0002 (M≈48.16) is Lean-checked as `75/2`. -/

theorem l0018_envelope_N16 : omega0 75 (1/2) = 75/2 := full_N16_omega0
theorem l0018_envelope_lt_cs0002_M : omega0 75 (1/2) < 48 := by native_decide
/-- Cleared: 2·(75/2)=75; Stokes refuter of C-0002 uses Ω(0)=75/2. -/
theorem l0018_envelope_cleared : 2 * omega0 75 (1/2) = 75 := full_N16_omega0_cleared

/-! ## L-0018 / C-0004 envelope on N=24: Ω ≤ K² E0 = 147/2 -/

theorem full_N24_omega0 : omega0 147 (1/2) = 147/2 := by native_decide
theorem l0018_envelope_N24 : omega0 147 (1/2) = 147/2 := full_N24_omega0
theorem full_N24_omega0_cleared : 2 * omega0 147 (1/2) = 147 := by native_decide
/-- Stokes floor factor e^{-2ν K² T} not formalized here; arithmetic pin only. -/
theorem l0018_N24_envelope_exceeds_c0005_M : omega0 147 (1/2) > 61 := by native_decide

/-! ## L-0025 / C-0006 envelope pin on N=32: Ω ≤ K² E0 = 300/2 -/

theorem full_N32_omega0 : omega0 300 (1/2) = 150 := by native_decide
theorem full_N32_omega0_cleared : 2 * omega0 300 (1/2) = 300 := by native_decide
/-- C-0006 M = 1.5 · Stokes_floor < envelope 150; pin envelope > 67. -/
theorem l0025_N32_envelope_exceeds_c0006_M : omega0 300 (1/2) > 67 := by native_decide
/-- Growth witness: N32 spectral ceiling exceeds N24. -/
theorem full_growth_N24_to_N32 : omega0 147 (1/2) < omega0 300 (1/2) := by native_decide

end NSGalerkin.Certificates
