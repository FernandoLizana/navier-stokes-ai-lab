/-!
# L-0003 comparison-ODE algebraic closure (finite Nat)

Formalizes the *algebraic* steps that close the L-0003 uniform bound, without
claiming a full ODE existence theory in Lean:

1. Sign of the comparison force via product inequalities:
   - if `a ≤ b·z` then `a·(z·z) ≤ (b·z)·(z·z)` (dissipation dominates),
   - if `b·z ≤ a` then `(b·z)·(z·z) ≤ a·(z·z)` (stretch dominates).
2. Algebraic identity for `Ω_eq · (ν·ν) = 6·M·(E0·E0)`.
3. Uniform cap: `max(Ω0, Ω_eq)` dominates both.

**Honesty.** Nat algebra only. No continuum NS. No Clay.
Evidence: N8 skeleton for comparison-closure arithmetic.
-/

namespace NSGalerkin.ComparisonODE

def dissipDominates (a b z : Nat) : Prop := a ≤ b * z
def stretchDominates (a b z : Nat) : Prop := b * z ≤ a

/-- Dissipation dominates ⇒ `a (z²) ≤ (b z) (z²)`. -/
theorem force_nonpos_of_dissip
    (a b z : Nat) (h : dissipDominates a b z) :
    a * (z * z) ≤ (b * z) * (z * z) :=
  Nat.mul_le_mul_right (z * z) h

/-- Stretch dominates ⇒ `(b z) (z²) ≤ a (z²)`. -/
theorem force_nonneg_of_stretch
    (a b z : Nat) (h : stretchDominates a b z) :
    (b * z) * (z * z) ≤ a * (z * z) :=
  Nat.mul_le_mul_right (z * z) h

/-- Exact Nat equilibrium `a = b · z_eq` ⇒ both domination predicates. -/
theorem equilibrium_of_exact_div (a b z_eq : Nat) (h : a = b * z_eq) :
    dissipDominates a b z_eq ∧ stretchDominates a b z_eq := by
  constructor
  · show a ≤ b * z_eq
    exact h ▸ Nat.le_refl (b * z_eq)
  · show b * z_eq ≤ a
    exact h ▸ Nat.le_refl (b * z_eq)

/-- Hypothesis form of `Ω_eq · ν² = 6 M E0²`. -/
theorem omega_eq_product_form (M E0 nu Omega_eq : Nat)
    (h : Omega_eq * (nu * nu) = 6 * M * (E0 * E0)) :
    Omega_eq * (nu * nu) = 6 * M * (E0 * E0) := h

theorem a_sq_is_six_M (a M : Nat) (h : a * a = 6 * M) : a * a = 6 * M := h

theorem uniform_cap_ge_Omega0 (Omega0 Omega_eq : Nat) :
    Omega0 ≤ max Omega0 Omega_eq :=
  Nat.le_max_left Omega0 Omega_eq

theorem uniform_cap_ge_Omega_eq (Omega0 Omega_eq : Nat) :
    Omega_eq ≤ max Omega0 Omega_eq :=
  Nat.le_max_right Omega0 Omega_eq

/--
Composite structural statement used by CERT-L0003:
`Ω_cap = max(Ω0, Ω_eq)` with the product identity for `Ω_eq`.
-/
theorem uniform_cap_closes
    (Omega0 Omega_eq M E0 nu : Nat)
    (hΩ : Omega_eq * (nu * nu) = 6 * M * (E0 * E0)) :
    Omega0 ≤ max Omega0 Omega_eq
      ∧ Omega_eq ≤ max Omega0 Omega_eq
      ∧ Omega_eq * (nu * nu) = 6 * M * (E0 * E0) :=
  ⟨uniform_cap_ge_Omega0 Omega0 Omega_eq,
   uniform_cap_ge_Omega_eq Omega0 Omega_eq,
   omega_eq_product_form M E0 nu Omega_eq hΩ⟩

/-- Above equilibrium (`b z ≥ a`), the comparison force product is non-positive
    in the sense `a(z²) ≤ (bz)(z²)`. Together with the hand ODE comparison
    principle (N7, external), this is the algebraic half of uniform closure. -/
theorem above_eq_force_nonpos
    (a b z : Nat) (h : dissipDominates a b z) :
    a * (z * z) ≤ (b * z) * (z * z) :=
  force_nonpos_of_dissip a b z h

/-- Below or at equilibrium, stretch product dominates. -/
theorem below_eq_force_nonneg
    (a b z : Nat) (h : stretchDominates a b z) :
    (b * z) * (z * z) ≤ a * (z * z) :=
  force_nonneg_of_stretch a b z h

example : dissipDominates 3 1 3 ∧ stretchDominates 3 1 3 :=
  equilibrium_of_exact_div 3 1 3 rfl

example : 24 * (1 * 1) = 6 * 4 * (1 * 1) := by decide

end NSGalerkin.ComparisonODE
