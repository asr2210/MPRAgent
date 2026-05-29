# 005 — cCRE + motif-embedded synthetic (75k / 75k)

## Goal
Test whether motif content is what makes random useful as a diversity
source. Replaces the uniform-random half of exp 004 with random backgrounds
that have 2 JASPAR consensus motifs embedded.

## Method
- 75,000 cCRE-centered (seed=4)
- 75,000 synthetic: random ACGT background + 2 motifs embedded at random
  positions from JASPAR 2024 vertebrate non-redundant (~879 motifs;
  consensus = argmax base of PFM per position)
- Concatenated, shuffled, written.

## Result
Mean across 14 evals: **0.8829**. vs uniform-random hybrid (0.8804, +0.003).
vs pure cCRE (0.8862, −0.003). Runtime 3418s.

| eval | cCRE   | hybrid_uniform | hybrid_motif | Δ(motif−uniform) |
|------|--------|----------------|--------------|------------------|
| 01   | 0.8288 | 0.8217         | 0.8258       | +0.004 |
| 02   | 0.9270 | 0.9199         | 0.9232       | +0.003 |
| 03   | 0.9193 | 0.9116         | 0.9148       | +0.003 |
| 04   | 0.8657 | 0.8609         | 0.8608       | +0.000 |
| 05   | 0.8285 | 0.8214         | 0.8256       | +0.004 |
| 06   | 0.9274 | 0.9204         | 0.9235       | +0.003 |
| 07   | 0.9025 | 0.8896         | 0.8926       | +0.003 |
| 08   | 0.8922 | 0.9156         | 0.9160       | +0.000 |
| 09   | 0.9465 | 0.9387         | 0.9363       | −0.002 |
| 10   | 0.9272 | 0.9193         | 0.9242       | +0.005 |
| 11   | 0.8145 | 0.8078         | 0.8114       | +0.004 |
| 12   | 0.7959 | 0.7900         | 0.7926       | +0.003 |
| 13   | 0.9038 | 0.8885         | 0.8898       | +0.001 |
| 14   | 0.9276 | 0.9205         | 0.9237       | +0.003 |

## Key observations
1. **Motif-embedded is slightly better than uniform random as the diversity
   source** — +0.003 mean, consistent across most evals.
2. **Eval_08 is essentially the same** (0.9160 vs 0.9156). The diversity
   premium on eval_08 is NOT specifically about motif content — it comes
   from "the library contains non-cCRE sequences", regardless of whether
   those non-cCRE sequences are random ACGT or motif-rich. Important
   negative result.
3. **Still below pure cCRE on mean** (−0.003). 50/50 dilution is too much
   even with smarter diversity content.
4. **Hard evals (11, 12) tiny improvements** from motifs (+0.004, +0.003)
   — still 0.81 / 0.79.

## Theory update (v4 → v5)
- The eval_08 super-additivity is about **source diversity**, not motif
  content. Any non-cCRE sequence type (uniform random, motif-embedded
  random, presumably also dinucleotide-shuffled cCRE) unlocks that eval.
- Motif content gives a small but real bonus on most other evals (+0.003).
  Likely because motif-embedded sequences activate the model's "regulatory
  recognition" learning pathway more cleanly than uniform random does.
- Pure cCRE (0.886) is still the leader on mean. Hybrids ≤ 0.886. So
  there's no "free lunch" — diversity helps eval_08 but costs elsewhere.
- Hard evals (11, 12) are not unlocked by ANY diversity strategy tried so
  far. Probably need cell-type-specific elements, conservation-filtered
  cCREs, or activity-extreme sequences.

## Next
Two refinements to try:
- **A. 90/10 cCRE/motif-embedded** — refine the ratio. A small diversity
  component might preserve eval_08 win without paying full dilution cost,
  potentially beating 0.886.
- **B. cCRE class-balanced** — equal numbers of PLS, pELS, dELS-type and
  CTCF-type cCREs. Tests whether internal cCRE diversity matters.

Going with **A** — directly tests whether hybrids can ever beat pure cCRE
on mean. Cheap mechanistic test. If 90/10 wins, the design philosophy is
"mostly cCRE + small diversity". If 90/10 loses on mean (likely, given
trend), pure-cCRE is the ceiling for this family of strategies and we
need a new lever.
