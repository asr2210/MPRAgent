# 004 — cCRE + random hybrid (75k / 75k)

## Goal
Test the diversity hypothesis. Random had a unique strength on eval_08 that
cCREs only partially recovered. Does a 50/50 mix capture both?

## Method
- 75,000 cCRE-centered windows (same pool as exp 003), seed=3
- 75,000 uniform random 200bp ACGT sequences, seed=3
- Concatenated, shuffled, written.

## Result
Mean across 14 evals: **0.8804** (vs 0.8862 for pure cCRE, **−0.006**).
Runtime 4749s (slower — possibly contention).

| eval | random  | genome  | cCRE   | hybrid | Δ(hybrid−cCRE) |
|------|---------|---------|--------|--------|----------------|
| 01   | 0.7760  | 0.8182  | 0.8288 | 0.8217 | −0.007 |
| 02   | 0.8718  | 0.9162  | 0.9270 | 0.9199 | −0.007 |
| 03   | 0.8562  | 0.9080  | 0.9193 | 0.9116 | −0.008 |
| 04   | 0.8178  | 0.8567  | 0.8657 | 0.8609 | −0.005 |
| 05   | 0.7756  | 0.8181  | 0.8285 | 0.8214 | −0.007 |
| 06   | 0.8722  | 0.9166  | 0.9274 | 0.9204 | −0.007 |
| 07   | 0.7910  | 0.8882  | 0.9025 | 0.8896 | −0.013 |
| 08   | 0.9080  | 0.8679  | 0.8922 | **0.9156** | **+0.023** |
| 09   | 0.8890  | 0.9364  | 0.9465 | 0.9387 | −0.008 |
| 10   | 0.8670  | 0.9062  | 0.9272 | 0.9193 | −0.008 |
| 11   | 0.7623  | 0.8044  | 0.8145 | 0.8078 | −0.007 |
| 12   | 0.7394  | 0.7857  | 0.7959 | 0.7900 | −0.006 |
| 13   | 0.7762  | 0.8947  | 0.9038 | 0.8885 | −0.015 |
| 14   | 0.8722  | 0.9164  | 0.9276 | 0.9205 | −0.007 |

## Key observations
1. **Eval_08 jumps to 0.916** — higher than EITHER pure random (0.908) OR
   pure cCRE (0.892). The mix is super-additive on this eval. This is the
   diversity premium I was hoping to find, and it's real.
2. **Mean drops by 0.006** because of small losses on every other eval
   (−0.005 to −0.015). Diluting cCRE costs information density.
3. **Eval_07 and eval_13 lose most** (−0.013, −0.015). These are the evals
   where genome → cCRE gain was largest, so they're presumably the most
   "genomic context-dependent" — diluting half of those with random hurts
   most.
4. **Hard evals (11, 12) lose ≈ −0.006** — random doesn't help them either.
5. **Eval_08 is special.** Random is best, hybrid is even better, genome is
   worst. This eval rewards having both "functional sequences" AND
   "diverse / random-like sequences" in training. Maybe it tests prediction
   on near-baseline (low-activity, random-like) sequences where pure cCRE
   has under-represented the low-activity end of the distribution.

## Theory update (v3 → v4)
- Diversity is real but selective: it's super-additive on one eval and
  subtractive on others. The eval_08 effect is *not* a simple weighted
  average — it exceeds both single sources.
- This is consistent with: the MPRA produces a wide activity distribution
  on random; eval_08 may test models on sequences that are at the low /
  near-baseline activity range, which cCREs underrepresent.
- The right library probably has a *small* random-like component to anchor
  the model on baseline / low-activity, paired with a *large* cCRE
  component for high-information regulatory examples.
- Diluting cCRE 50/50 is too much. The right ratio is probably 80–90%
  cCRE + 10–20% random-like.

## Next
Two strong candidates:
- **A. 90/10 cCRE/random** — refines the hybrid ratio. Direct test of "small
  random component preserves eval_08 win without sacrificing cCRE density".
- **B. cCRE + motif-embedded synthetic** (75/75) — tests whether motif-rich
  synthetic outperforms uniform random as the diversity source. New lever.

Going with **B** — different axis, more conceptual info. If motif-embedded
beats uniform random as the diversity source, we've isolated motifs as the
active ingredient. If equal, random's contribution is about coverage / GC
diversity rather than motif content. Either result reshapes future designs.
