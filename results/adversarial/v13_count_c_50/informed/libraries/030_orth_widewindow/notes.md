# Experiment 030 — wider orth-DHS window (350bp)

## Result
eval_01 = 0.5753 (**−0.0009 vs 025's 0.5762**). Wider exclusion window
slightly hurt.

Pool shrank 904k → 604k (still 24x N_DHS, not coverage-limited).
Tightening orthogonality past the 200bp default lost more in DHS
diversity than it gained in orthogonality.

## Numbers
| eval | 025 (W=200) | **030 (W=350)** | Δ      |
|------|------------:|----------------:|-------:|
| 01   | 0.5762      | 0.5753          | −0.0009|
| 04/09| 0.5620      | 0.5564          | −0.006 |
| 07   | 0.6131      | 0.6175          | +0.004 |
| 08   | 0.1591      | 0.1419          | −0.017 |
| 13   | 0.5925      | 0.5980          | +0.005 |

Wider window pushes the pool further from cCRE-adjacent sequences,
which boosts the regulatory-grammar tasks (07/13) but costs the
cell-type-specific tasks (04/09) and the CpG-tinged eval_08.

## Conclusion
**025 (W=200) is the optimal orth-DHS window.** Per-eval tradeoffs net
out to −0.0009 on eval_01. The 200bp window precisely separates
"summit inside a cCRE" from "summit independently called as a DHS in
a quiet region" — exactly what the orthogonal pool needs.

## Final answer
After 30 experiments, the best library design is **025**:
- 25,000 orth-DHS sequences (Meuleman Index, nsamp ≤ 5, ≥200bp from
  nearest cCRE midpoint)
- 25,000 cCRE sequences balanced across 6 functionally-specific
  classes (CA-CTCF, CA-H3K4me3, CA-TF, PLS, TF, pELS — i.e.
  excluding dELS and CA)

eval_01 = 0.5762, mean across 14 evals ≈ 0.555.
