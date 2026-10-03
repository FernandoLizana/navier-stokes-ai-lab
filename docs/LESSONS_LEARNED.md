# Lessons learned (documented failures)

Sources: `reports/C0008_CERT_REPAIR_STATUS.md`, `reports/C0008_REPAIR_EVALUATION.md`, `reports/L0073R_CLUSTER_AUDIT.md`, and related repair code under `python/ns_exploration/terminal_weighted/` and `python/ns_exploration/validation/`.  
Each case is short: goal → defect → effect → reported fix → remaining limit.

---

### 1. Frobenius streaming without entry aggregation

| | |
|--|--|
| **Goal** | Upper-bound operator norms via Frobenius-type accumulation |
| **Defect** | Summing `|delta|^2` per update instead of aggregating into matrix entries then summing `|M_ij|^2` |
| **Effect** | Incorrect / non-comparable bounds |
| **Reported fix** | Entry aggregation then square (`C0008_CERT_REPAIR_STATUS.md` P0 #1) |
| **Evidence now** | Repair docs mark FIXED; full re-audit of every historical Frobenius artifact not repeated here |
| **Limit** | Float / MPFR paths still need their stated evidence tags |

### 2. Invalid cluster `max(φ)·B_C` inequality

| | |
|--|--|
| **Goal** | Tighter Route A' integral via clusters |
| **Defect** | Multiplying residual by outer `max_phi` (invalid for matrix sums); ~50× false deflation |
| **Effect** | Historical ~37.92 looked “better” but was not a valid repair bound |
| **Reported fix** | L-0073R residual reweight; audit shows outer factor removed (`L0073R_CLUSTER_AUDIT.md`) |
| **Evidence now** | Corrected cluster sum ≈ 38.52 reported; still does not close C-0008 |
| **Limit** | Triangle inequality across clusters remains an N7 gap |

### 3. Verifiers trusting declared totals / missing manifests

| | |
|--|--|
| **Goal** | Adversarial certificate verification |
| **Defect** | Trusting declared integrals; soft-fail on missing manifests |
| **Effect** | PASS could disagree with recomputation |
| **Reported fix** | Recompute from blocks; fail closed (`C0008_CERT_REPAIR_STATUS.md`) |
| **Limit** | PASS = checks implemented, not continuum truth |

### 4. Interval constants and sqrt / Omega orientation

| | |
|--|--|
| **Goal** | Rigorous-ish interval enclosures |
| **Defect** | Float constants; wrong sqrt orientation for Omega reported |
| **Effect** | Interval claims weakened or wrong |
| **Reported fix** | Documented FIXED in repair status |
| **Limit** | Exact coefficient enclosure still listed as pending in places |

### 5. Direct-vs-tensor test compared a quantity to itself

| | |
|--|--|
| **Goal** | Consistency between construction paths |
| **Defect** | Test compared `f_from_G(G,z)` to itself |
| **Effect** | False confidence |
| **Reported fix** | Four-path comparison (repair status P0 #5) |
| **Limit** | Always ask what independent quantity a test compares |

### 6. Permutation multiplicities collapsed by sets

| | |
|--|--|
| **Goal** | Streaming assembly of interaction terms |
| **Defect** | Using sets of permutations drops multiplicities |
| **Effect** | Wrong combinatorial weights |
| **Reported fix** | Explicit 6-tuple lists (repair status P0 #6) |
| **Limit** | Combinatorial bugs are silent without enumeration tests |

### 7. Historical “best” vs repair honesty vs closure

| | |
|--|--|
| **Goal** | Close C-0008 (I_hi < I_* ≈ 1.30) |
| **Defect** | Treating historical cluster ~37.92 as rigorous; conflating verifier PASS with mathematical closure |
| **Effect** | Misleading narrative of near-success |
| **Reported fix** | Primary Route A ≈ 170.59 (honest, far); cluster sketch ≈ 38.52; status remains `exploring` |
| **Limit** | Still ~30–131× above the sharp integral target under current methods |

---

**Takeaway:** software PASS, interval tags, and Lean stubs are different layers. Documented bugs show why aggressive optimism is dangerous in AI-assisted math tooling.
