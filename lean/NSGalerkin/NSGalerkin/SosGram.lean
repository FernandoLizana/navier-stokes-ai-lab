/-!
# N5 SOS Gram certificate arithmetic (Lean N8 bridge, L-0059)

Cross-checks the conservative rational upper bounds `C_ub_hi` from
CERT-L0052/53/54 against `C_dagger` using `native_decide`.

**Honesty.** Arithmetic of exported certificate constants only. Does NOT
formalize TSSOS, Gram PSD, or the cubic stretch inequality. No Clay implication.
-/

namespace NSGalerkin.SosGram

/-- Low-slab stretch target C_dagger (conservative rational upper bound). -/
def C_dagger : Rat := 9562009938690147 / 10^15

/-! ## CERT-L0052 band {1,2,3} D=52 -/

def l0052_C_ub_hi : Rat := 843822489307524 / 10^15

theorem l0052_C_ub_lt_Cdagger : l0052_C_ub_hi < C_dagger := by native_decide

/-! ## CERT-L0054 band {1,2,3,4} D=64 -/

def l0054_C_ub_hi : Rat := 894040132411397 / 10^15

theorem l0054_C_ub_lt_Cdagger : l0054_C_ub_hi < C_dagger := by native_decide

/-! ## CERT-L0053 one-pol C-R-0008 D=92 -/

def l0053_C_ub_hi : Rat := 3356065194637164 / 10^15

theorem l0053_C_ub_lt_Cdagger : l0053_C_ub_hi < C_dagger := by native_decide

/-! ## CERT-L0058 band {1,2,3,4,5} D=112 -/

def l0058_C_ub_hi : Rat := 1305261518419503 / 10^15

theorem l0058_C_ub_lt_Cdagger : l0058_C_ub_hi < C_dagger := by native_decide

/-! ## CERT-L0060 band {1,2,3,4,5,6} D=112 (interval constant) -/

def l0060_C_ub_hi : Rat := 1700932221482412 / 10^15

theorem l0060_C_ub_lt_Cdagger : l0060_C_ub_hi < C_dagger := by native_decide

/-! ## CERT-L0061 band {1..8} D=184 (COSMO interval) -/

def l0061_C_ub_hi : Rat := 1851834363071286 / 10^15

theorem l0061_C_ub_lt_Cdagger : l0061_C_ub_hi < C_dagger := by native_decide

/-! ## L-0062 COSMO ladder {1..9}, {1..10}, {1..12} -/

def l0062_band9_C_ub_hi : Rat := 2223426008475950 / 10^15
def l0062_band10_C_ub_hi : Rat := 2440195403421871 / 10^15
def l0062_band12_C_ub_hi : Rat := 2782805118332654 / 10^15

theorem l0062_band9_C_ub_lt_Cdagger : l0062_band9_C_ub_hi < C_dagger := by native_decide
theorem l0062_band10_C_ub_lt_Cdagger : l0062_band10_C_ub_hi < C_dagger := by native_decide
theorem l0062_band12_C_ub_lt_Cdagger : l0062_band12_C_ub_hi < C_dagger := by native_decide

/-! ## CERT-L0065 one-pol greedy-second D=154 -/

def l0065_C_ub_hi : Rat := 4014938261573984 / 10^15

theorem l0065_C_ub_lt_Cdagger : l0065_C_ub_hi < C_dagger := by native_decide

/-! ## CERT-L0067 one-pol greedy-second+24 D=178 -/

def l0067_C_ub_hi : Rat := 40760825002657395 / 10^16

theorem l0067_C_ub_lt_Cdagger : l0067_C_ub_hi < C_dagger := by native_decide

/-- Eleven structured N5 constant certs sit strictly below C_dagger. -/
theorem all_eleven_below_Cdagger :
    l0052_C_ub_hi < C_dagger ∧
    l0054_C_ub_hi < C_dagger ∧
    l0053_C_ub_hi < C_dagger ∧
    l0058_C_ub_hi < C_dagger ∧
    l0060_C_ub_hi < C_dagger ∧
    l0061_C_ub_hi < C_dagger ∧
    l0062_band9_C_ub_hi < C_dagger ∧
    l0062_band10_C_ub_hi < C_dagger ∧
    l0062_band12_C_ub_hi < C_dagger ∧
    l0065_C_ub_hi < C_dagger ∧
    l0067_C_ub_hi < C_dagger := by
  constructor
  · exact l0052_C_ub_lt_Cdagger
  constructor
  · exact l0054_C_ub_lt_Cdagger
  constructor
  · exact l0053_C_ub_lt_Cdagger
  constructor
  · exact l0058_C_ub_lt_Cdagger
  constructor
  · exact l0060_C_ub_lt_Cdagger
  constructor
  · exact l0061_C_ub_lt_Cdagger
  constructor
  · exact l0062_band9_C_ub_lt_Cdagger
  constructor
  · exact l0062_band10_C_ub_lt_Cdagger
  constructor
  · exact l0062_band12_C_ub_lt_Cdagger
  constructor
  · exact l0065_C_ub_lt_Cdagger
  · exact l0067_C_ub_lt_Cdagger

/-- Legacy alias (ten certs before L-0067). -/
theorem all_ten_below_Cdagger :
    l0052_C_ub_hi < C_dagger ∧
    l0054_C_ub_hi < C_dagger ∧
    l0053_C_ub_hi < C_dagger ∧
    l0058_C_ub_hi < C_dagger ∧
    l0060_C_ub_hi < C_dagger ∧
    l0061_C_ub_hi < C_dagger ∧
    l0062_band9_C_ub_hi < C_dagger ∧
    l0062_band10_C_ub_hi < C_dagger ∧
    l0062_band12_C_ub_hi < C_dagger ∧
    l0065_C_ub_hi < C_dagger := by
  constructor
  · exact l0052_C_ub_lt_Cdagger
  constructor
  · exact l0054_C_ub_lt_Cdagger
  constructor
  · exact l0053_C_ub_lt_Cdagger
  constructor
  · exact l0058_C_ub_lt_Cdagger
  constructor
  · exact l0060_C_ub_lt_Cdagger
  constructor
  · exact l0061_C_ub_lt_Cdagger
  constructor
  · exact l0062_band9_C_ub_lt_Cdagger
  constructor
  · exact l0062_band10_C_ub_lt_Cdagger
  constructor
  · exact l0062_band12_C_ub_lt_Cdagger
  · exact l0065_C_ub_lt_Cdagger

end NSGalerkin.SosGram
