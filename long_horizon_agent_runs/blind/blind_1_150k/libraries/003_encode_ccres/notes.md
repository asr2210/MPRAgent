# 003 — ENCODE cCREs (200bp centered)

## Goal
Test the function-enrichment axis: do sequences centered on putative
regulatory elements (ENCODE V3 cCREs) give further gains beyond plain
genomic tiling?

## Method
- ENCODE V3 GRCh38 cCREs (1,063,878 elements across all classes and cell types)
- For each cCRE: 200bp window centered on (start+end)/2
- Same filter as exp 002: skip if OOB, any non-ACGT, >50% softmasked
- 745,398 windows survived; 150,000 sampled uniformly (seed=2)
- cCRE class distribution unmodified (natural mix: ~48% dELS, ~26% dELS+CTCF, ~9% pELS+CTCF, etc.)

## Result vs prior experiments
| eval | random  | genome  | cCRE   | Δ(cCRE−genome) |
|------|---------|---------|--------|----------------|
| 01   | 0.7760  | 0.8182  | 0.8288 | +0.0106 |
| 02   | 0.8718  | 0.9162  | 0.9270 | +0.0108 |
| 03   | 0.8562  | 0.9080  | 0.9193 | +0.0113 |
| 04   | 0.8178  | 0.8567  | 0.8657 | +0.0090 |
| 05   | 0.7756  | 0.8181  | 0.8285 | +0.0104 |
| 06   | 0.8722  | 0.9166  | 0.9274 | +0.0108 |
| 07   | 0.7910  | 0.8882  | 0.9025 | +0.0143 |
| 08   | 0.9080  | 0.8679  | 0.8922 | **+0.0243** |
| 09   | 0.8890  | 0.9364  | 0.9465 | +0.0101 |
| 10   | 0.8670  | 0.9062  | 0.9272 | +0.0210 |
| 11   | 0.7623  | 0.8044  | 0.8145 | +0.0101 |
| 12   | 0.7394  | 0.7857  | 0.7959 | +0.0102 |
| 13   | 0.7762  | 0.8947  | 0.9038 | +0.0091 |
| 14   | 0.8722  | 0.9164  | 0.9276 | +0.0112 |

Aggregate: mean **0.8862**. K562 0.8841, HepG2 0.8933, SKNSH 0.8812.
Time 3041s (faster than 002).

## Key observations
1. **cCRE beats genome on every eval, including the eval_08 regression.**
   No losses anywhere — function enrichment is a strict win.
2. **Eval_08 jumped most** (+0.0243) — the regression in exp 002 is largely
   recovered. This is interesting: eval_08 specifically penalises plain
   genomic context but rewards focused regulatory content. Random (0.908) is
   still slightly better on eval_08 than cCRE (0.892), so eval_08 must be
   testing something where breadth of activity range / wide sequence
   coverage helps.
3. **Diminishing returns on function enrichment.** random→genome was +0.053
   on mean, genome→cCRE only +0.013. The function-axis is saturating.
4. **Hard evals still hard.** eval_11 (0.815) and eval_12 (0.796) are the
   floor, only +0.01 over genome. Whatever they test isn't unlocked by
   focusing on regulatory regions.
5. **HepG2 gained most** across all three experiments (0.820 → 0.877 → 0.893).
   K562 and SKNSH track closely. No cell-line systematically left behind.

## Theory update (v2 → v3)
- The function-enrichment axis explains a lot — ~+0.07 from random to cCRE.
  But it's plateauing.
- "Diversity helps" hypothesis (from eval_08 in exp 002) is partially
  supported: cCREs are STRICTLY better than genome (no losses), but they
  still don't fully recover random's eval_08 strength. Something in random
  is still useful that cCREs lack.
- The hard evals (11, 12) suggest there's a class of test sequences not
  well-represented in any of: random, genome, cCRE. Maybe they test
  cell-type-specific biology, extreme activity values, or non-cCRE
  regulatory elements (e.g., insulators in non-CTCF contexts, polycomb
  domains, or evolutionarily young regulatory elements not annotated by
  ENCODE).

## Next
Two strongest candidates for experiment 004:
- **A. cCRE + random hybrid** (75k each). Tests the "diversity helps"
  hypothesis directly. If mean > 0.886 *and* eval_08 > 0.892, diversity
  beats pure-function and the next several experiments should explore
  multi-source libraries.
- **B. Synthetic library with embedded JASPAR motifs in random
  background.** Tests "explicit motifs in controlled context" — different
  axis. Could be very informative if it competes with cCREs.

Going with **A** — it's the direct test of the eval_08 signal we already see,
and it informs the basic design philosophy (pure-function vs diverse) which
shapes everything downstream.
