import NSGalerkin.EnstrophyEnergy

/-!
# Discrete enstrophy–energy inequality over nonnegative rationals

Lifts the Nat lemma `enstrophy_le_K2_energy` to mode energies in `ℚ≥0`
represented with a **common positive denominator** `D`:

  eᵢ = nᵢ / D,   nᵢ ∈ ℕ,   D ≠ 0.

Multiplying through by `D` reduces the claim to the Nat inequality already
proved. This matches Fourier mode energies `|û_k|² ≥ 0` written with a common
clearing denominator (w.l.o.g. for any finite list of rationals).

**Honesty.** Finite lists only; no continuum NS; no Clay. Evidence: N8 skeleton.
-/

namespace NSGalerkin.EnstrophyEnergyRat

open NSGalerkin.EnstrophyEnergy

abbrev ModeNums := List (Nat × Nat)

/-- Rational energy of a mode: `n / D`. -/
def energyRat (D n : Nat) : Rat := (n : Rat) / (D : Rat)

/-- Weighted enstrophy sum over ℚ: `Σ k² · (n/D)`. -/
def weightSumRat (D : Nat) : ModeNums → Rat
  | [] => 0
  | (k2, n) :: rest => (k2 : Rat) * energyRat D n + weightSumRat D rest

/-- Energy sum over ℚ: `Σ (n/D)`. -/
def energySumRat (D : Nat) : ModeNums → Rat
  | [] => 0
  | (_k2, n) :: rest => energyRat D n + energySumRat D rest

theorem mul_energyRat (D n : Nat) (hD : (D : Rat) ≠ 0) :
    (D : Rat) * energyRat D n = (n : Rat) := by
  -- D * (n / D) = n
  simpa [energyRat, Rat.mul_comm] using Rat.div_mul_cancel (a := (n : Rat)) hD

theorem mul_weightSumRat (D : Nat) (hD : (D : Rat) ≠ 0) (xs : ModeNums) :
    (D : Rat) * weightSumRat D xs = (weightSum xs : Rat) := by
  induction xs with
  | nil =>
      simp [weightSumRat, weightSum]
  | cons p rest ih =>
      rcases p with ⟨k2, n⟩
      -- D * (↑k2 * energyRat + rest) = ↑k2 * (D * energyRat) + D * rest
      simp only [weightSumRat, weightSum]
      rw [Rat.mul_add]
      -- first term: D * (↑k2 * energyRat) = ↑k2 * (D * energyRat) = ↑k2 * ↑n
      have h1 :
          (D : Rat) * ((k2 : Rat) * energyRat D n)
            = (k2 : Rat) * (n : Rat) := by
        calc
          (D : Rat) * ((k2 : Rat) * energyRat D n)
              = (k2 : Rat) * ((D : Rat) * energyRat D n) := by
                rw [← Rat.mul_assoc, Rat.mul_comm (D : Rat) (k2 : Rat), Rat.mul_assoc]
          _ = (k2 : Rat) * (n : Rat) := by rw [mul_energyRat D n hD]
      rw [h1, ih]
      -- ↑k2 * ↑n + ↑(weightSum rest) = ↑(k2*n + weightSum rest)
      rw [← Rat.natCast_mul k2 n, ← Rat.natCast_add (k2 * n) (weightSum rest)]

theorem mul_energySumRat (D : Nat) (hD : (D : Rat) ≠ 0) (xs : ModeNums) :
    (D : Rat) * energySumRat D xs = (energySum xs : Rat) := by
  induction xs with
  | nil =>
      simp [energySumRat, energySum]
  | cons p rest ih =>
      rcases p with ⟨k2, n⟩
      simp only [energySumRat, energySum]
      rw [Rat.mul_add, mul_energyRat D n hD, ih]
      rw [← Rat.natCast_add n (energySum rest)]

/--
Main ℚ lemma: nonnegative mode energies with common den `D ≠ 0`,
`Σ k²·(n/D) ≤ K² · Σ(n/D)` whenever every `k² ≤ K²`.
-/
theorem enstrophy_le_K2_energy_rat
    (K2 D : Nat) (hD : (D : Rat) ≠ 0) (hDpos : 0 < (D : Rat))
    (xs : ModeNums) (h : boundedBy K2 xs) :
    weightSumRat D xs ≤ (K2 : Rat) * energySumRat D xs := by
  have hNat := enstrophy_le_K2_energy K2 xs h
  have hW := mul_weightSumRat D hD xs
  have hE := mul_energySumRat D hD xs
  have hCast : (weightSum xs : Rat) ≤ (K2 : Rat) * (energySum xs : Rat) := by
    have h1 : (weightSum xs : Rat) ≤ ((K2 * energySum xs : Nat) : Rat) :=
      (Rat.natCast_le_natCast).mpr hNat
    calc
      (weightSum xs : Rat) ≤ ((K2 * energySum xs : Nat) : Rat) := h1
      _ = (K2 : Rat) * (energySum xs : Rat) := Rat.natCast_mul K2 (energySum xs)
  have hScaled :
      (D : Rat) * weightSumRat D xs
        ≤ (D : Rat) * ((K2 : Rat) * energySumRat D xs) := by
    calc
      (D : Rat) * weightSumRat D xs
          = (weightSum xs : Rat) := hW
      _ ≤ (K2 : Rat) * (energySum xs : Rat) := hCast
      _ = (K2 : Rat) * ((D : Rat) * energySumRat D xs) := by rw [hE]
      _ = (D : Rat) * ((K2 : Rat) * energySumRat D xs) := by
          rw [← Rat.mul_assoc, Rat.mul_comm (K2 : Rat) (D : Rat), Rat.mul_assoc]
  exact Rat.le_of_mul_le_mul_left hScaled hDpos

/-- Helper: `D ≠ 0` as Nat ⇒ `↑D ≠ 0` and `0 < ↑D`. -/
theorem nat_cast_pos_of_ne_zero (D : Nat) (h : D ≠ 0) :
    (D : Rat) ≠ 0 ∧ 0 < (D : Rat) :=
  ⟨fun h0 => h (Rat.natCast_eq_zero_iff.mp h0), Rat.natCast_pos.mpr (Nat.pos_of_ne_zero h)⟩

/-- Convenient form from a Nat `D ≠ 0` hypothesis. -/
theorem enstrophy_le_K2_energy_rat'
    (K2 D : Nat) (hD : D ≠ 0) (xs : ModeNums) (h : boundedBy K2 xs) :
    weightSumRat D xs ≤ (K2 : Rat) * energySumRat D xs := by
  have ⟨hne, hpos⟩ := nat_cast_pos_of_ne_zero D hD
  exact enstrophy_le_K2_energy_rat K2 D hne hpos xs h

/-- Concrete smoke: K²=4, D=2, numerators giving energies 1/2, 3/2, 1. -/
example :
    weightSumRat 2 [(1, 1), (4, 3), (0, 2)]
      ≤ (4 : Rat) * energySumRat 2 [(1, 1), (4, 3), (0, 2)] := by
  apply enstrophy_le_K2_energy_rat'
  · decide
  · intro p hp
    simp at hp
    rcases hp with rfl | rfl | rfl <;> decide

end NSGalerkin.EnstrophyEnergyRat
