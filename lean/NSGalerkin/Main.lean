import NSGalerkin

/-!
Executable smoke entry: confirms the library builds.
Does not claim continuum NS or Clay results.
-/

def main : IO Unit := do
  IO.println "NSGalerkin OK:"
  IO.println "  - EnstrophyEnergy: discrete Ω ≤ K²·E (Nat)"
  IO.println "  - EnstrophyEnergyRat: same over ℚ≥0 with common denominator"
  IO.println "  - ComparisonODE: L-0003 force sign + uniform cap algebra"
  IO.println "  - Certificates: shell-k3 / N12 / N16 exact caps + growth"
  IO.println "Scope: finite-mode identities only. NOT continuum NS. NOT Clay."
