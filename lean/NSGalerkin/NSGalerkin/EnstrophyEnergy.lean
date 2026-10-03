/-!
# Discrete enstrophy–energy inequality (finite modes)

Formalizes the elementary combinatorial bound used as the first step of
L-0001…L-0006:

  If each retained mode satisfies `k² ≤ K²` and `e ≥ 0` is a mode energy
  weight, then
      Σ (k² · e) ≤ K² · Σ e.

In the continuum notation of the project this is `Ω ≤ K² E` (the factor ½
cancels on both sides).

**Scope (mathematical honesty):**
* This is a statement about **finite lists of natural numbers**.
* It does **not** formalize the Navier–Stokes PDE, Leray projection, or Clay.
* Evidence target: N8 *skeleton* for one elementary identity only.
-/

namespace NSGalerkin.EnstrophyEnergy

/-- Weighted enstrophy sum: `Σ k² · e`. -/
def weightSum : List (Nat × Nat) → Nat
  | [] => 0
  | (k2, e) :: rest => k2 * e + weightSum rest

/-- Energy sum: `Σ e`. -/
def energySum : List (Nat × Nat) → Nat
  | [] => 0
  | (_k2, e) :: rest => e + energySum rest

/-- Every mode's squared wavenumber is at most `K2`. -/
def boundedBy (K2 : Nat) (xs : List (Nat × Nat)) : Prop :=
  ∀ p, p ∈ xs → p.1 ≤ K2

/--
Main lemma (discrete `Ω ≤ K² E`):
for any finite list of mode pairs `(k², e)` with `k² ≤ K2`,
`weightSum xs ≤ K2 * energySum xs`.
-/
theorem enstrophy_le_K2_energy
    (K2 : Nat) (xs : List (Nat × Nat)) (h : boundedBy K2 xs) :
    weightSum xs ≤ K2 * energySum xs := by
  induction xs with
  | nil =>
      simp [weightSum, energySum]
  | cons p rest ih =>
      have hp : p.1 ≤ K2 := h p (by simp)
      have hrest : boundedBy K2 rest := fun q hq => h q (List.mem_cons_of_mem p hq)
      have ih' : weightSum rest ≤ K2 * energySum rest := ih hrest
      have hterm : p.1 * p.2 ≤ K2 * p.2 := Nat.mul_le_mul_right p.2 hp
      calc
        weightSum (p :: rest)
            = p.1 * p.2 + weightSum rest := rfl
        _ ≤ K2 * p.2 + weightSum rest := Nat.add_le_add_right hterm _
        _ ≤ K2 * p.2 + K2 * energySum rest := Nat.add_le_add_left ih' _
        _ = K2 * (p.2 + energySum rest) := by rw [Nat.mul_add]
        _ = K2 * energySum (p :: rest) := rfl

/-- Empty support is trivial. -/
theorem enstrophy_le_K2_energy_nil (K2 : Nat) :
    weightSum [] ≤ K2 * energySum [] := by
  simp [weightSum, energySum]

/-- Concrete smoke check: modes with K² = 4. -/
example : weightSum [(1, 1), (4, 2), (0, 3)] ≤ 4 * energySum [(1, 1), (4, 2), (0, 3)] := by
  apply enstrophy_le_K2_energy
  intro p hp
  simp at hp
  rcases hp with rfl | rfl | rfl <;> decide

end NSGalerkin.EnstrophyEnergy
