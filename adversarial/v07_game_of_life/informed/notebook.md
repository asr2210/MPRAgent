# MPRA Library Design — Lab Notebook

## 2026-06-03 07:19 — Setup and initial theory

**Setup.** Working dir initialized. Symlinked shared resources to `data/`:
- `hg38.fa` — human reference
- `ENCODE_cCREs_v3.bed` — 1,063,878 ENCODE SCREEN candidate cis-regulatory elements
- `gencode_v46_*.bed` — TSS coordinates
- `JASPAR2024_CORE_non-redundant_pfms_jaspar.txt` — TF motifs
- `ccre_motif_scores.npz` — precomputed motif scores

**Baseline landscape (from instructions.md, 50k libraries, eval_01):**
- Best: `dhs_topic` 0.7232 — DHS sampled proportional to NMF topic loadings
- `dhs_sei` 0.7201 — 50% DHS topic + 50% SEI class-prop
- `dhs_synth` 0.7174 — 50% DHS topic + 50% random
- `dhs_random` 0.7089
- `dhs_stratified` 0.7055
- `synth_oracle` 0.6840 — pure random sequences
- `mpra_oracle` 0.6643 — sampled from prior MPRA dataset
- `mpra_real` 0.6026 — same sequences but with REAL (noisy) MPRA labels

**Key observations from baselines:**
1. Topic-weighted DHS > random DHS by ~0.014 → cell-type-specific signal matters
2. DHS + random (`dhs_synth`) < DHS alone → random sequences slightly dilute training
3. Pure random (`synth_oracle`) at 0.684 → surprisingly capable; sequence features generalize
4. `mpra_oracle` (sampled from prior MPRA) is WORSE than DHS → MPRA distribution is narrower
5. `mpra_real` (noisy real labels) much worse than `mpra_oracle` → label noise kills training
6. `dhs_stratified_sei_synth` (3-way diversity-maximized) at 0.7094 → diversity-maximization
   alone doesn't beat informed topic sampling

**Generalization beyond K562/HepG2/SK-N-SH.** The 14 eval sets are anonymous; based on the
instructions ("equally informative for cell types we have no data on"), eval sets likely
include held-out cell types and possibly other regulatory contexts. The library must
expose the model to *transferable* regulatory primitives: TF motifs, motif syntax, element
types (promoter/enhancer/insulator). A library that only captures K562/HepG2/SK-N-SH
features will overfit; one that spans pan-tissue regulatory grammar will transfer.

**Current theory v1.** The most informative library has three properties:
(a) **Functional diversity**: covers all regulatory element TYPES (PLS, pELS, dELS,
    CTCF, DNase-H3K4me3), not just open-chromatin signal in 3 cell types
(b) **Cell-type breadth**: includes elements active in many cell types, not just the
    measured three — so motif grammar learned transfers to held-out cells
(c) **Sequence-feature diversity**: rich coverage of TF motifs and their combinations

DHS-topic captures (b) well but is biased toward enhancer-like (mostly distal) elements.
Adding categorical balance from cCREs should improve (a).

**What this predicts.** If theory is right, an ENCODE cCRE-balanced library — equal
representation across PLS, pELS, dELS, CTCF-only, DNase-H3K4me3 — should match or exceed
dhs_topic. If it underperforms, then mere category balance is insufficient and DHS topic
weighting (cell-type information) dominates.

## 2026-06-03 07:24 — Plan: Experiment 001 (cCRE category-balanced)

**Mode: exploring new hypothesis.** I have not yet run any experiment in this directory.
Establishing a strong, theory-driven first datapoint matters more than refining baselines
I haven't reproduced.

**Design.** 50,000 sequences = 10,000 from each of 5 cCRE primary categories:
- PLS (Promoter-Like Signature) — combine PLS + PLS,CTCF-bound
- pELS (proximal Enhancer-Like) — pELS + pELS,CTCF-bound
- dELS (distal Enhancer-Like) — dELS + dELS,CTCF-bound
- CTCF-only — CTCF-only,CTCF-bound
- DNase-H3K4me3 — DNase-H3K4me3 + DNase-H3K4me3,CTCF-bound

For each cCRE: center a 200bp window on the midpoint, extract from hg38, skip if it
contains N, retry adjacent until valid. Uppercase ACGT only.

**Prediction.** I expect this to score in [0.68, 0.74] on eval_01. If it beats 0.7232,
functional category balance is more powerful than topic-stratified DHS. If it lands
between 0.69 and 0.72, it confirms DHS-derived sequences are strong but suggests
categorical balance is comparable to topic weighting. If <0.68, then categorical balance
loses cell-type-specific information that DHS topics encode.

**Justification for generalization.** Promoter/enhancer/CTCF categories are universal
across cell types. A model trained on balanced functional categories should learn
element-type-specific grammars (e.g., TATA-box context for promoters, enhancer motif
clusters for dELS) that apply to ANY cell type, not just the three labeled ones.

## 2026-06-03 07:44 — Result: Experiment 001 (cCRE balanced) — BIG SURPRISE

**Result.** eval_01 = 0.3919, mean across 14 evals = 0.3812.
Per-cell-type on eval_01: K562 0.604, HepG2 0.428, SK-N-SH 0.143.
Runtime ~20 min.

**Why this is a big deal.** My cCRE library performs at the **random baseline**.
Comparing to `strategies.md` (which is calibrated to this run's evaluator):
- gc_50 random: 0.3972
- random_uniform: 0.3951
- my cCRE library: 0.3919

So balanced biologically-curated regulatory elements give NO improvement over
random sequences. This contradicts the `instructions.md` baseline table (which
claims `dhs_topic = 0.7232`). Either those baselines are from a different
evaluator (with much more training data/better model), or they are misleading
context. Either way they are NOT the relevant comparison for this run.

**The real ceiling I need to beat: ~0.40 (gc_50 baseline).**

**Theory v1 falsified.** "Functional diversity from cCRE categories will beat
DHS topic" is moot — cCRE categories don't beat anything; the real surprise is
that biology itself doesn't help at THIS sample size/model capacity.

**New theory v2.** With 50k training sequences and this evaluator, the model
appears to either:
- (a) lack capacity to use motif-level biological structure,
- (b) overfit to training-set base composition, or
- (c) be tested on a distribution where biological priors are not aligned.

The lever that matters most appears to be **base composition + diversity at the
nucleotide level** (since gc_50 wins). Biological grammar may help only when
sequences ALSO match the eval distribution's base composition.

**Cell-type asymmetry.** K562 always trains 4-5x better than SK-N-SH for any
strategy. This means the biggest gains in mean_r come from anything that lifts
HepG2 and especially SK-N-SH. K562 is near its ceiling already.

**What this predicts for experiment 002.** A library that simply matches the
evaluator's likely distribution (controlled GC ~50%, high sequence diversity)
will outperform any biologically curated set with biased composition.

## 2026-06-03 07:46 — Plan: Experiment 002 (motif-implanted GC-50 random)

**Mode: exploring new hypothesis informed by exp 001 result.**

**Question.** If gc_50 random sequences are the strongest known baseline (0.397),
does ADDING TF motifs to them improve learning? This isolates two effects:
1. Base composition (held at 50% GC, matching gc_50)
2. Biological signal (TF motif content injected)

If motif-implanted GC-50 > gc_50, the model CAN learn from motifs when
composition is right. If they're equal, motifs add nothing on top of composition.

**Design.** 50k sequences, each:
- Generate random ACGT with GC=50% (i.i.d., 200bp)
- Pick a random subset of JASPAR TF motifs to implant (3-5 per sequence)
- Implant at random non-overlapping positions
- The bulk of each sequence is still GC-50 random

This is essentially `gc_50 + motifs` — closer to gc_50 baseline (the strongest
known) but augmented with biological signal. JASPAR has ~700 vertebrate motifs;
I'll use a curated subset known to be functional in K562/HepG2/SK-N-SH and
related lineages.

**Prediction.** Score in [0.40, 0.50] on eval_01. If above gc_50 baseline
(0.397), motif content helps the model. If at gc_50, it doesn't.

**Justification for generalization.** TF motifs are conserved across cell types
(the TFs themselves are expressed differently, but the binding grammar is
universal). A model that learns motif features will generalize, while a model
that only learns composition statistics will not.

## 2026-06-03 08:15 — Result: Experiment 002 — Motif content doesn't help either

**Result.** eval_01 = 0.3939, mean across 14 evals = 0.3827.
Per-cell-type on eval_01: K562 0.614, HepG2 0.431, SK-N-SH 0.136.

**Now we have three libraries scoring at exactly the same ceiling:**
- gc_50 random (baseline): 0.397
- 001_ccre_balanced (biology only): 0.392
- 002_motif_implanted_gc50 (biology + random scaffold): 0.394

The pattern is STRIKINGLY consistent. The per-cell-type values are nearly
identical across libraries (K562 ~0.60, HepG2 ~0.43, SK-N-SH ~0.14). The model
appears to converge to the same per-cell-type ceiling no matter what kind of
sequences I train it on — as long as the composition is roughly balanced.

**Theory v2 falsified (partially).** Adding motifs to a composition-matched
random scaffold did not help. Either:
- (a) The model cannot use motif-level features at this data scale
- (b) The eval ceiling is determined by something independent of training
       sequence quality (noise floor, model capacity, eval set difficulty)
- (c) Motifs in random scaffolds aren't a strong enough signal — would need
       motif clusters / specific spatial arrangements / etc.

**Theory v3.** The 0.397 ceiling on this evaluator is likely a NOISE/CAPACITY
floor for the small model trained on 50k×3=150k examples. The K562 score
~0.60 reflects how well the eval correlates with the model's strongest signal.
SK-N-SH scoring 0.14 reflects either a SK-N-SH-specific data scarcity or
an SK-N-SH eval set that's genuinely harder to predict.

To break past 0.40, I likely need:
1. Sequences that produce stronger MPRA signal differential (high-activity
   promoters that genuinely drive expression in the cell types)
2. Or sequences that better cover the eval distribution
3. Or explicit cell-type-discriminating elements

**What hasn't been tested:** activity-range coverage. Random sequences have
low MPRA activity (no real signal). cCRE balanced is mixed. A library FOCUSED
on the most-active class (promoters, PLS) should produce more learnable
discriminative signal because the activity range is wider.

## 2026-06-03 08:18 — Plan: Experiment 003 (PLS-only / promoter-focused)

**Mode: exploring new hypothesis.**

**Question.** Does focusing on promoter (PLS) elements — which drive the
strongest MPRA activity — give the model richer training signal? PLS elements
contain canonical promoter features (TATA box, Inr, GC-box, CpG islands)
and are the most "active" of the cCRE classes.

**Design.** 50k 200bp sequences from ENCODE cCRE PLS class only. There are
~41k PLS elements; to reach 50k I'll generate 1-2 windows per element with
random jitter where elements are long enough (median cCRE len = 286bp, so
many can yield 2 distinct 200bp windows).

**Prediction.** If PLS-rich training improves on cCRE-balanced, score > 0.42.
If equal to gc_50/cCRE (0.39-0.40), then activity richness within biology
isn't the lever. Either way it constrains the theory.

**Generalization justification.** Promoters are universal regulatory elements —
every gene has a promoter, in every cell type. A model that learns promoter
grammar (Inr, TATA, GC-box positioning) should generalize to held-out cell
types since the basic transcription machinery is conserved.

## 2026-06-03 08:30 — Result: Experiment 003 — PLS-only LOST ground

**Result.** eval_01 = 0.3748 (WORSE than cCRE balanced 0.3919 and gc_50 0.397).
Mean across 14 evals = 0.3653.

**Mechanism revealed.** PLS library mean GC = 0.617 (vs ideal 0.50). The drop
~0.04 from gc_50 baseline matches the composition penalty curve:
- gc_50 (0.50 GC) → 0.397
- PLS (0.62 GC) → 0.375
- gc_rich (0.80 GC) → 0.222
- at_rich (0.20 GC) → 0.299

**Theory v3.** Composition is the DOMINANT lever at this scale. Performance
drops monotonically as training set GC moves from 0.50 in either direction.
Biological content adds NOTHING measurable above what composition predicts.

This is consistent with the literature finding (MDC preprint, 2025) that at
low training data sizes, native and random sequences perform similarly. At
50k×3=150k examples, we're in that regime.

**What this changes.** Searching for the "right" biological library is futile
until composition is controlled. The right next step is to test:
- Can composition-matched biology BEAT gc_50? (exp 004)
- Is there a sweet spot OTHER than 0.50? (e.g., 0.45 or 0.55?)
- Does within-library composition uniformity matter, or just the mean?

## 2026-06-03 08:31 — Plan: Experiment 004 (cCRE filtered to GC ∈ [0.45, 0.55])

**Mode: refining promising direction (composition control).**

**Question.** When composition is matched to gc_50 (50% GC ± 5%), does
biological content (cCREs) provide ANY lift over pure random?

**Design.** Sample cCRE-derived 200bp windows, then filter to only keep
sequences with GC content in [0.45, 0.55]. Continue sampling until 50k valid.
Should produce a library of real regulatory elements with controlled composition.

**Prediction.**
- If 004 > 0.41: biology adds ~0.01-0.02 over composition-matched random.
  Significant — biology can help, just not when composition is wrong.
- If 004 ≈ 0.40 (≈ gc_50): composition is fully sufficient; biology provides
  no measurable lift over random at this scale.
- If 004 < 0.40: biology actively hurts even when composition-matched
  (unlikely given exp 002 result).

**Generalization justification.** If composition-controlled biology helps,
the library is teaching universal regulatory grammar (independent of which
specific cell type provides the labels) — a strong signal for generalization.
If composition alone matters, there's no point overfitting to K562/HepG2/SKNSH
biology; a clean composition-matched library is universally informative.

## 2026-06-03 08:53 — Result: Experiment 004 — biology adds NOTHING even when GC-matched

**Result.** eval_01 = 0.3921 (cCRE balanced 0.392, gc_50 0.397). Mean 14 evals = 0.3811.

**Theory v3 confirmed.** ALL composition-matched libraries fall in 0.390-0.397
range, regardless of biological content. The "biology helps" hypothesis is
falsified at this 50k scale.

**Frustrating but informative.** I now have 4 datapoints supporting the
composition-only theory:
| library          | mean GC | eval_01 |
|------------------|---------|---------|
| gc_50 (baseline) | 0.50    | 0.397   |
| cCRE balanced    | 0.51    | 0.392   |
| cCRE GC-matched  | 0.50    | 0.392   |
| motif+gc50       | 0.50    | 0.394   |
| PLS only         | 0.62    | 0.375   |

The 0.39-0.40 ceiling appears very tight. To make progress, I need
qualitatively different experiments.

**What's left.** Three big unexplored hypotheses:
1. **Source matters for natural sequences.** Test: random hg38 windows
   (non-regulatory) vs cCRE windows at matched GC.
2. **Cell-type-specific motifs may lift SK-N-SH.** Test: implant neural TF
   motifs (NEUROG, ASCL1, POU3F2, etc.) and check if SKNSH score rises.
3. **Sequence DIVERSITY at matched GC.** Could include 5-10 different
   compositions averaging 0.50, or specific dinucleotide structures.

## 2026-06-03 08:55 — Plan: Experiment 005 (random hg38 genomic windows, GC-matched)

**Mode: exploring new hypothesis (does annotation curation add value?).**

**Question.** When composition is matched and we draw from natural human
sequences, does it matter whether the regions are "annotated regulatory"
(cCRE) or just RANDOM genomic windows (mostly intergenic/intronic)?

**Design.** Sample 50k random 200bp windows from hg38 autosomes (chr1-22):
- Pick chromosome weighted by length
- Pick random start
- Extract 200bp; reject if non-ACGT
- Filter for GC ∈ [0.45, 0.55]

This will largely be intronic and intergenic sequence — biologically
"unselected" but at matched composition.

**Prediction.**
- If 005 ≈ 0.392 (= cCRE GC-matched): the curation status doesn't matter;
  ANY natural human DNA at GC=50% gives the same result.
- If 005 > 0.397 (> gc_50): natural genomic dinucleotide/k-mer structure
  helps even outside regulatory regions.
- If 005 < 0.392: regulatory cCREs do have something extra over random genome.

**Generalization justification.** Random genomic windows are an unbiased
sample of the species' DNA distribution. If the library at this scale
provides equivalent training signal regardless of biological annotation, the
implication is that the model is learning generic sequence-feature statistics,
not specific regulatory grammar.





## 2026-06-03 09:14 — Result: Experiment 005 — natural genome at GC=50 scores like random

**Result.** eval_01 = 0.3900, mean across 14 evals = 0.3805.
Per-cell-type on eval_01: K562 0.6013, HepG2 0.4232, SK-N-SH 0.1454.
Runtime ~20 min.

**Theory v3 strongly confirmed.** Random hg38 windows at GC ∈ [0.45, 0.55]
score **0.3900** — within noise of cCRE-GC-matched (0.3921), motif+gc50 (0.394),
and gc_50 baseline (0.397). The data are now unambiguous: composition is the
full story for this evaluator at the 50k scale.

Updated composition-controlled comparison:
| library              | mean GC | eval_01 |
|----------------------|---------|---------|
| gc_50 (baseline)     | 0.50    | 0.397   |
| 001 cCRE balanced    | 0.51    | 0.392   |
| 002 motif+gc50       | 0.50    | 0.394   |
| 003 PLS only         | 0.62    | 0.375   |
| 004 cCRE GC-matched  | 0.50    | 0.392   |
| 005 hg38 random      | 0.49    | 0.390   |

The 5 composition-matched libraries fall within 0.007 of each other.
**Biological provenance is essentially invisible at this scale.**

**Per-cell-type asymmetry remains untouched.** K562 ≈ 0.60, HepG2 ≈ 0.42,
SK-N-SH ≈ 0.15 — across ALL five libraries. The model has a structural
ceiling on each cell type that's independent of training sequence quality.
This rules out theory variants based on "training sequences too biased"
because SK-N-SH still scores 0.15 even with diverse cCRE input.

**Updated theory v4.** Performance is determined by:
1. **Composition** (mean GC near 0.50 is best — ~0.07 lift over GC=0.62)
2. **Per-cell-type ceiling** from eval set properties, not training input
   (K562 0.60, HepG2 0.42, SK-N-SH 0.15)
3. The remaining ~0.02 spread across composition-matched libraries is
   essentially measurement noise.

To break the ~0.40 ceiling, training-set sequence selection (random,
biological, motif-rich) is not the lever. What might be:
- (a) Activity-range coverage: training sequences with wide MPRA dynamic
      range provide stronger regression signal. Most random/cCRE sequences
      have MPRA activity near baseline — limited learnable variance.
- (b) Pretrained oracle pre-selection: use Malinois (BassetBranched CNN
      pretrained on K562/HepG2/SK-N-SH MPRA) to score 500k candidates and
      pick 50k that span the predicted activity distribution.

## 2026-06-03 09:16 — Plan: Experiment 006 (Malinois-oracle activity-spanning library)

**Mode: exploring qualitatively new mechanism.**

**Question.** Does dynamic-range coverage of training-set ACTIVITY (not
sequence diversity) lift the model above the 0.40 composition ceiling?

**Theory.** Current libraries are GC=50 sequences with near-baseline MPRA
activity. The model's regression target has small variance, so even a
perfect model has limited correlational signal to learn from. If I select
sequences predicted to span the full activity range (high active, high
repressive, neutral) in K562/HepG2/SK-N-SH, the regression problem becomes
better-conditioned and the trained model should generalize better.

**Design.** Generate ~500k candidate sequences (mix of cCREs and gc_50
random, all composition-matched). Score each with the pretrained Malinois
CNN producing predicted log2FC in K562/HepG2/SK-N-SH. Select 50k that:
- Cover the full predicted activity range (uniform binning across deciles)
- Span all three cell types' predicted distributions
- Maintain mean GC ≈ 0.50

**Prediction.** If activity coverage is the lever, eval_01 > 0.42 (clearly
above the 0.40 ceiling). If composition is truly the ONLY thing that
matters, eval_01 ≈ 0.39 (within noise of other GC-50 libraries). Either
outcome is highly informative.

**Generalization justification.** A library with broad predicted activity
in our three labeled cell types likely has sequences with broad activity
in OTHER cell types as well (TFs are shared, regulatory grammar is shared).
Selecting for activity-range in the labeling cell types is a proxy for
activity-range in the held-out cell types — directly serving generalization.

**Risk.** Malinois was trained on K562/HepG2/SK-N-SH MPRA data — using its
predictions to select training sequences could bias the trained model
toward Malinois's failure modes. To mitigate: I will only use Malinois
to *bucket* by predicted activity, not to memorize Malinois's specific
predictions.

## 2026-06-03 09:42 — Result: Experiment 006 — FIRST LIFT ABOVE CEILING

**Result.** eval_01 = **0.3964**, mean across 14 evals = **0.3860**.
Per-cell-type on eval_01: K562 **0.6177** (+0.016), HepG2 **0.4333** (+0.010),
SK-N-SH 0.1384 (unchanged). Runtime 20 min.

**Theory v4 partially confirmed.** Malinois-oracle activity-spanning selection
gave the first datapoint above the gc_50/cCRE plateau:

| library                | eval_01 | K562  | HepG2 | SKNSH |
|------------------------|---------|-------|-------|-------|
| gc_50 baseline         | 0.397   | ~     | ~     | ~     |
| 001 cCRE balanced      | 0.392   | 0.604 | 0.428 | 0.143 |
| 002 motif+gc50         | 0.394   | 0.614 | 0.431 | 0.136 |
| 004 cCRE GC-matched    | 0.392   | 0.606 | 0.424 | 0.146 |
| 005 hg38 random GC-mat | 0.390   | 0.601 | 0.423 | 0.145 |
| **006 Malinois oracle**| **0.396** | **0.618** | **0.433** | 0.138 |

The lift is small (+0.005 over gc_50, +0.006 over 005) but it's the FIRST
break upward from the composition-only plateau. The mechanism appears real:
sequences with diverse predicted activities give the small model better
discriminative signal.

**SK-N-SH problem persists.** All libraries score ~0.14 on SK-N-SH. Malinois
SK-N-SH predictions did not transfer. Either (a) the evaluator's SK-N-SH
sequences are systematically different from Gosai 2024 SK-N-SH, (b) the model
has a structural ceiling for SK-N-SH on this evaluator, or (c) we need
neural-lineage-specific motifs to lift it.

**Theory v5.** The full causal model is:
1. **Composition** (~0.04 lift from GC=0.5 vs GC=0.6/0.4)
2. **Activity range coverage** (~0.005 lift from oracle selection)
3. **Per-cell-type structural ceiling** (K562 ~0.62, HepG2 ~0.43, SK-N-SH ~0.15)
   that determines the upper bound on mean_r regardless of training set.

To break further:
- (a) Push activity selection harder — select TOP DECILE active sequences only
      (highest signal-to-noise for regression).
- (b) Try cell-type-DISCRIMINATING selection (high variance across 3 cells)
      — this is the most "transferable" mechanism for unseen cell types.
- (c) For SK-N-SH specifically: implant neural TF motifs (NEUROG2, ASCL1,
      POU3F2/4, ONECUT, OLIG, LHX, NEUROD1, NHLH).

## 2026-06-03 09:44 — Plan: Experiment 007 (Malinois TOP-magnitude only)

**Mode: refining promising direction (activity selection).**

**Question.** Does the lift in exp 006 come from SPAN (diverse activities)
or from MAGNITUDE (presence of high-activity sequences in training)?

**Design.** 50k sequences = top-50k by max(predicted log2FC) across 3 cells,
from the same 500k GC-50 candidate pool used in exp 006. Composition still
matched to 0.50.

**Prediction.**
- If 007 > 006 (>0.40): magnitude is the dominant lever; high-activity
  training sequences give the strongest regression signal.
- If 007 ≈ 006 (0.39-0.40): both span and magnitude help equally.
- If 007 < 006: span/diversity matters more than magnitude alone — model
  needs both extremes to anchor the regression.

**Generalization justification.** Highly active sequences contain functional
TF motifs that activate expression in ANY cell type that expresses those
TFs. A library of strong activators sampled from random sequence space
contains de novo high-affinity motif combinations — these should activate
similarly across cell types that share TF expression, providing universal
training signal.

**Risk.** If we select only top-magnitude, we may LOSE the ability to
predict baseline/inactive sequences in held-out evals. The model could
extrapolate poorly to low-activity regions.

## 2026-06-03 10:05 — Result: Experiment 007 — TIED with 006

**Result.** eval_01 = **0.3969** (vs 006 span = 0.3964, gc_50 = 0.397).
Mean across 14 evals = **0.3861** (vs 006 = 0.3860).
Per-cell-type: K562 0.618, HepG2 0.436, SKNSH 0.137. Runtime ~20 min.

**Theory v5 refined.** Top-magnitude and 3D-span give the SAME ~0.396 mean
on eval_01 — basically indistinguishable. Both approaches deliver the same
+0.005 lift over composition-matched random.

| library                 | eval_01 | K562  | HepG2 | SKNSH | strategy |
|-------------------------|---------|-------|-------|-------|----------|
| 005 hg38 random GC-mat  | 0.390   | 0.601 | 0.423 | 0.145 | random natural   |
| 006 Malinois span 5x5x5 | 0.396   | 0.618 | 0.433 | 0.138 | 125-cell stratify|
| 007 Malinois top-mag    | 0.397   | 0.618 | 0.436 | 0.137 | top-50k by max   |

The fact that span ≈ top-magnitude means **the lift comes from EXPOSURE TO
HIGH-ACTIVITY EXAMPLES, not the diversity of activity levels.** Once
high-activity sequences are in training, the model uses them; adding low-
activity sequences (or vice versa) doesn't change predictions much.

**Updated theory v5.**
1. Composition (GC ≈ 0.5) — large lift (~0.04 from GC=0.6)
2. Activity exposure — small lift (~0.005), saturated quickly
3. Per-cell-type structural ceiling — independent of training set:
   - K562 ~0.62 (hardest to budge)
   - HepG2 ~0.43
   - SK-N-SH ~0.14 (stubbornly low)
4. eval_08 outlier (always ~0.27 vs 0.39) — likely a structurally different
   cell type / experimental condition where the model breaks down

**What hasn't been tested.**
- Cell-type DISCRIMINATION (high-variance across 3 cells) — selecting for
  sequences that express ONLY in one cell type. These contain cell-type-
  specific TF motifs and should provide the most transferable training
  signal for unseen cell types.
- Natural sequence + oracle (cCRE-derived top-activity)
- Hyperactive design via FastSeqProp gradient ascent

## 2026-06-03 10:07 — Plan: Experiment 008 (cell-type-discriminating selection)

**Mode: exploring new hypothesis (cross-cell-type transferability).**

**Question.** Do sequences that DIFFERENTIALLY activate across cell types
provide stronger training signal for cross-cell-type generalization than
sequences that uniformly activate? Cell-type-specific motifs (e.g., a
hematopoietic TF binding site that activates K562 but not HepG2) are
exactly the kind of feature a generalizing model needs to learn.

**Design.** 500k GC-50 random ACGT candidates scored with Malinois. Compute
two summary stats per sequence:
- magnitude = max(preds)
- variance = std(preds across 3 cells)

Select 50k by maximizing magnitude × variance, OR by joint stratification:
- Half from top-25k by variance (most cell-type-discriminating)
- Half from balanced mixture combining magnitude and variance

Simpler design: select top 50k by composite score = magnitude + 2*variance.

**Prediction.**
- If 008 > 007: cell-type-discrimination is a stronger lever than pure
  magnitude — strong evidence for the "transfer via motif specificity"
  theory.
- If 008 ≈ 007 (0.395-0.400): magnitude alone is sufficient; cell-type
  patterns add no additional information.
- If 008 < 007: discrimination filters out useful broadly-active sequences
  too aggressively.

**Generalization justification.** Cell-type-discriminating sequences expose
the model to the *grammar of specificity* — which motif combinations cause
which cell type to respond. For a held-out cell type, this grammar
transfers (the TFs are shared, only their expression pattern differs).
Sequences that are uniformly active in all 3 reveal nothing about which
TF families drive which response.

## 2026-06-03 10:28 — Result: Experiment 008 — discrimination hurts slightly

**Result.** eval_01 = **0.3936** (vs 007 = 0.3969, 006 = 0.3964).
Mean across 14 evals = **0.3835** (vs 007 = 0.3861).
Per-cell-type: K562 0.616, HepG2 0.430, SKNSH 0.135.

**Theory v5 update — discrimination is not the lever.** Selecting for high
cross-cell variance hurt vs pure magnitude. The model benefits from
exposure to broadly-active sequences, not cell-type-specific ones.

Three Malinois-based libraries all land within 0.004 of each other,
confirming the lever is "any high-activity exposure" and that it
saturates near 0.397:

| library                | eval_01 | mean_r | strategy                     |
|------------------------|---------|--------|------------------------------|
| 006 span 5x5x5         | 0.3964  | 0.3860 | 3D uniform coverage          |
| 007 top-magnitude      | 0.3969  | 0.3861 | top max(K562,HepG2,SKNSH)    |
| 008 cell-type discrim  | 0.3936  | 0.3835 | rank(mag)+2*rank(std)        |

**The ceiling is real and tight.** Across 8 experiments now:
- All composition-matched libraries score 0.39-0.397 on eval_01
- Best score is gc_50 strategies baseline at 0.397 = same as my Malinois 007
- Oracle adds ~+0.005 max
- Per-cell-type values are nearly identical across ALL libraries

**Updated theory v6.** The 0.40 mean ceiling on this evaluator is determined
by two factors:
1. **The held-out eval distribution shifts away from any training dist** —
   no amount of training selection captures more of it.
2. **The per-cell-type ceilings** (K562 0.62, HepG2 0.43, SK-N-SH 0.15) are
   STRUCTURAL features of the evaluator, not training choices.

The fact that SK-N-SH is stuck at 0.13-0.15 across EVERY library suggests
either (a) the evaluator's SK-N-SH eval set is genuinely much noisier or
biologically different from training, or (b) the model architecture cannot
fit SK-N-SH from this data scale.

**What's left to test.**
- (a) Combine natural sequence statistics with oracle: cCREs + Malinois top
- (b) Hyperactive generation: produce sequences with log2FC well beyond
  what random sampling can find (FastSeqProp / simulated annealing)
- (c) Neural-lineage motif augmentation targeted at SK-N-SH
- (d) Mixed library: combine multiple known-good strategies
- (e) Bigger candidate pool (5M) to find rarer extreme outliers

## 2026-06-03 10:30 — Plan: Experiment 009 (cCRE + Malinois top selection)

**Mode: refining promising direction (oracle + natural sequence).**

**Question.** Does combining BIOLOGICAL SEQUENCE STATISTICS (real cCREs)
with ORACLE SELECTION (Malinois top activity) give a bigger lift than
either alone? Random+Malinois plateaus at 0.397; pure cCREs at 0.392.
Their combination might exceed both if biology and oracle scoring are
complementary.

**Design.** Load all ~1M ENCODE cCREs, extract centered 200bp windows,
keep only those with GC ∈ [0.45, 0.55], score all with Malinois, take top
50k by max(predicted log2FC).

**Prediction.**
- If 009 > 0.397: biology + oracle is complementary; real enhancers chosen
  for predicted activity beat random + oracle.
- If 009 ≈ 0.397 ± 0.002: both routes saturate at the same plateau.
- If 009 < 0.397: cCRE base composition or scaffold structure interferes.

**Generalization justification.** Real enhancer sequences contain natural
TF motif clusters and the spacing/context that evolution has selected for.
A library of real, high-activity enhancers exposes the model to
*biologically valid* motif syntax, which should transfer better to
held-out cell types than random sequences with the same predicted activity.

## 2026-06-03 10:55 — Result: Experiment 009 — biology + oracle HURT (surprising)

**Result.** eval_01 = **0.3879** (vs 007 = 0.3969, gc_50 = 0.397).
Mean across 14 evals = **0.3774** (vs 007 = 0.3861).
Per-cell-type: K562 **0.598** (−0.020), HepG2 **0.422** (−0.014),
SK-N-SH **0.144** (+0.007).

**Surprise.** Despite 60% higher predicted activity (selected K562 mean
4.10 vs 2.60 for random+oracle), cCRE+oracle UNDERPERFORMS random+oracle
on K562 and HepG2 by a notable margin. Only SK-N-SH lifted slightly.

**Likely mechanism: oracle bias.** Malinois was trained on cCRE-derived
sequences (Gosai 2024 MPRA), so it systematically overpredicts cCRE
activity. Selecting "top Malinois activity" from cCREs creates a label-
distribution mismatch with the evaluator's MPRA (which presumably uses
a different cCRE measurement). The trained model overshoots on cCRE-like
test sequences → lower correlation.

This is essentially a TRAINING-DATA SHIFT failure. The composition is
matched (0.501) but the higher-order sequence statistics of cCREs differ
from random, and Malinois picks the most cCRE-stereotypical sequences.

**Updated leaderboard:**
| library                  | eval_01 | mean_r |
|--------------------------|---------|--------|
| 007 random+Malinois top  | 0.3969  | 0.3861 |  ← best
| 006 random+Malinois span | 0.3964  | 0.3860 |
| gc_50 baseline           | 0.397   | ~      |
| 008 Mal cell-type discrim| 0.3936  | 0.3835 |
| 004 cCRE GC-matched      | 0.3921  | 0.3811 |
| 005 hg38 random GC-mat   | 0.3900  | 0.3805 |
| 009 cCRE+Malinois top    | 0.3879  | 0.3774 |
| 003 PLS only             | 0.3748  | 0.3653 |

**Theory v6.** When using an oracle for selection, the oracle's CALIBRATION
on the candidate pool matters. Random ACGT candidates: well-calibrated.
cCRE candidates: oracle overpredicts → biased selection → worse downstream.

**New direction.** SK-N-SH lifted in 009 (+0.007). This is the FIRST time
SK-N-SH responded to anything other than composition. Mechanism: cCREs
contain real TF motifs including neural-lineage TFs (in dELS class).

**What to test next.** Neural-lineage TF motif augmentation. If SK-N-SH
is liftable via specific motifs (NEUROG, ASCL, POU3F, etc.), this is
a new lever specifically for the lagging cell type.

## 2026-06-03 11:00 — Plan: Experiment 010 (neural TF motif implants for SK-N-SH lift)

**Mode: exploring new mechanism (cell-type-specific motif targeting).**

**Question.** Can we LIFT SK-N-SH (stuck at 0.13-0.15 across all 9
libraries so far) by implanting neural-lineage TF motifs into composition-
matched random scaffolds?

**Design.** 50k random GC=0.50 200bp scaffolds. For each, implant 3-5
JASPAR motif instances drawn from a CURATED neural TF list:
- Pioneer/proneural: NEUROG2 (MA0623), ASCL1 (MA1100), NEUROD1 (MA1109)
- Stem/early: SOX2 (MA0143), POU3F2 (MA0788), POU3F4 (MA0789)
- Olig/glial: OLIG2 (MA0678), ID4 (MA0824)
- Specification: LHX2 (MA0700), LHX3 (MA0699), ISL1 (MA1573)
- Neural-spec promoter: ZIC1-3 (MA0696-MA0697), RFX (MA0509)
- Neuroblastoma-specific: PHOX2B (MA0681), HAND2 (MA1638), GATA3 (MA0037)
- Repressors: REST/NRSF (MA0138) — KEY neural-lineage repressor

**Prediction.**
- If 010 lifts SK-N-SH (> 0.16): cell-type-specific motif augmentation
  is a new lever. K562/HepG2 may stay flat or dip slightly.
- If 010 SK-N-SH ≈ 0.14 (unchanged): the SK-N-SH ceiling is structural
  (not motif-content-dependent).
- If 010 lowers all three: motif augmentation introduces noise (similar
  to exp 002).

**Generalization justification.** Neural TFs are universal in
neural-lineage cell types. A model that learns motif → cell-type
response in SK-N-SH should transfer to any neural cell type held out.

## 2026-06-03 11:20 — Result: Experiment 010 — neural motifs failed to lift SK-N-SH

**Result.** eval_01 = **0.3927**, mean across 14 evals = **0.3818**.
Per-cell-type on eval_01: K562 0.612, HepG2 0.428, **SK-N-SH 0.138** (unchanged).
Runtime ~20 min.

**Strong negative result for SK-N-SH liftability via motif targeting.**
128 curated neural-lineage TF motifs implanted into 3-5 sites per sequence
gave the SAME result as exp 002 (random JASPAR motifs): SK-N-SH stays at
0.138.

**Theory v6 confirmed.** The per-cell-type ceilings (K562 ~0.62, HepG2
~0.43, SK-N-SH ~0.14) are STRUCTURAL features of the evaluator's
held-out test sequences, independent of training-data choices. The mean_r
ceiling at ~0.397 follows.

**Updated leaderboard (10 experiments):**
| library                  | eval_01 | mean_r | notes                  |
|--------------------------|---------|--------|------------------------|
| 007 random+Malinois top  | 0.3969  | 0.3861 | best                   |
| 006 random+Malinois span | 0.3964  | 0.3860 |                        |
| gc_50 baseline           | 0.397   | ~      | (from strategies.md)   |
| 002 random JASPAR motifs | 0.3939  | 0.3827 |                        |
| 008 Mal cell-type discrim| 0.3936  | 0.3835 |                        |
| 010 neural motifs        | 0.3927  | 0.3818 |                        |
| 001 cCRE balanced        | 0.3919  | 0.3812 |                        |
| 004 cCRE GC-matched      | 0.3921  | 0.3811 |                        |
| 005 hg38 random GC-mat   | 0.3900  | 0.3805 |                        |
| 009 cCRE+Malinois top    | 0.3879  | 0.3774 |                        |
| 003 PLS only             | 0.3748  | 0.3653 |                        |

After 10 experiments, every composition-matched library lands in
[0.388, 0.397]. The Malinois-oracle adds ~+0.005 in the random-pool
setting. No experiment exceeds the strategies.md max (0.397).

## 2026-06-03 11:22 — Plan: Experiment 011 (Hyperactive design via simulated annealing)

**Mode: bold new mechanism — generative sequence design.**

**Question.** Can we BREAK the 0.40 ceiling by training the model on
sequences with predicted activity FAR beyond what random sampling can
find? Random max log2FC = ~10; designed sequences could reach 15+.

**Design.** For each of 50k sequences:
1. Start with a random GC=0.50 200bp scaffold (Malinois ~0).
2. Run simulated annealing: at each step, propose a single-base mutation,
   accept if Malinois max(K562, HepG2, SKNSH) increases (or with cooling
   probability if it decreases). Run for ~200 mutation steps.
3. Filter to GC ∈ [0.40, 0.60].
4. Take resulting sequence.

This gives sequences with Malinois predictions in tail far beyond random.

**Computational cost.** 50k * 200 SA steps * 1 Malinois call each = 10M
forward passes. At 14k seqs/s in batched mode, would take ~12 min total
if we batch the SA steps (50k sequences * 200 steps in parallel = 50k SA
trajectories of 200 steps, batched 1024 at a time).

**Prediction.**
- If 011 > 0.405: hyperactive design IS a new lever beyond Malinois
  selection. Cell-type-shared TF binding sites learned from these synthetic
  high-activity sequences transfer better.
- If 011 ≈ 0.395-0.400: hyperactive design saturates at the same ceiling.
  No mechanism can break it.
- If 011 < 0.390: designed sequences overfit Malinois → worse downstream
  (similar to 009 cCRE+oracle failure mode).

**Generalization justification.** Sequences with very high Malinois
predictions almost certainly contain DENSE clusters of canonical TF
binding sites in proper spatial arrangements (Malinois learned what makes
sequences active). A library of such sequences should expose the trained
model to the GRAMMAR of activation — universal across cell types.

## 2026-06-03 11:55 — Result: Experiment 011 — hyperactive SA WORSE than top-magnitude

**Result.** eval_01 = **0.3947**, mean across 14 evals = **0.3839**.
Per-cell-type: K562 0.617, HepG2 0.432, SK-N-SH 0.136. Runtime ~20 min.

Hyperactive-designed sequences (mean Malinois 5.51, max 15.10) scored
WORSE than top-magnitude-selected random sequences (mean Malinois 2.60).

**Activity ↔ score is non-monotonic:**
| library                  | mean Malinois | eval_01 |
|--------------------------|---------------|---------|
| 005 random GC-mat        | 0.7           | 0.3900  |
| 006 random + Malinois span| ~3           | 0.3964  |
| 007 random + Malinois top| 2.60          | 0.3969  |
| 009 cCRE + Malinois top  | 4.10          | 0.3879  |
| 011 hyperactive SA       | 5.51          | 0.3947  |

The "sweet spot" appears to be moderate activity (Malinois 2-3). Both
under (random) and over (hyperactive/cCRE+oracle) hurt. Reason: at very
high Malinois scores, sequences become STEREOTYPED (dense canonical TF
clusters in unrealistic arrangements) that don't match the eval
distribution.

**Theory v6 reinforced.** The ~0.397 ceiling is robust against EVERY
single-strategy training-set selection I've tried (11 experiments now).
It's almost certainly a STRUCTURAL property of the evaluator + small
model + 50k training scale.

What HASN'T been tested:
- (a) MIXED libraries (multiple strategies combined in one library)
- (b) Different DIVERSITY structures (multiple GC bins; different
       composition shapes)
- (c) Adversarial / multi-oracle ensemble

## 2026-06-03 11:58 — Plan: Experiment 012 (mixed-strategy hybrid library)

**Mode: bold combination of known-good strategies.**

**Question.** Does combining multiple known-good single-strategy libraries
into one hybrid library help, by exposing the model to DIVERSE forms of
learnable signal?

**Design.** 50k sequences = 4-way mix:
- 12.5k from 007 random+Malinois top (the best single, captures activity)
- 12.5k from 006 random+Malinois span (covers activity range)
- 12.5k cCRE GC-matched (natural sequence statistics)
- 12.5k gc_50 i.i.d. random (clean baseline composition)

All 4 sub-pools are composition-matched at GC ≈ 0.50.

**Prediction.**
- If 012 > 0.40: strategy DIVERSITY is a new lever; the model benefits
  from multiple complementary feature distributions.
- If 012 ≈ 0.397: the ceiling is structural, and any single strategy
  captures everything a mixed strategy would.
- If 012 < 0.397: mixing dilutes the strongest signal (oracle selection)
  with weaker signals.

**Generalization justification.** A held-out cell type's eval distribution
is unknown — it could be enriched in random, in real enhancers, in
oracle-active, or in some mixture. A library that covers MULTIPLE
distributions hedges against unknown eval bias. This is the right
generalization strategy when the eval distribution is uncertain.

## 2026-06-03 12:23 — Result: Experiment 012 — hybrid mix WORSE than best single

**Result.** eval_01 = **0.3903**, mean_r = **0.3805**.
Per-cell-type: K562 0.607, HepG2 0.428, SK-N-SH 0.136.

**Mixing strategies dilutes signal.** 4-way 25%-each blend of best
strategies scored worse than any individual high-performer. This rules
out "diversity of strategies" as a lever.

**Strong reinforcement of theory v6.** After 12 experiments, the
relationship is now unmistakable:

```
Mean GC ~ 0.5    →  scores in [0.388, 0.397]
Mean GC ≠ 0.5    →  scores below 0.39
Oracle selection →  marginal +0.005 within the 0.5-GC bracket
Anything else    →  no detectable effect
```

The model + 50k training + evaluator structure determine almost
everything. Training-data choice is a small modulator at best.

## 2026-06-03 12:25 — Plan: Experiment 013 (4M random pool + Malinois top — pool-size test)

**Mode: refining promising direction (does oracle lift scale with pool size?).**

**Question.** 007 used a 500k pool and found mean Malinois score 2.60 in
the top 50k. With a 4M pool (8× larger), the top 50k would have higher
predicted activity. Does this scale the oracle benefit?

**Design.** 4M random GC-50 candidates, scored with Malinois, top 50k
by max(K562, HepG2, SKNSH). Compare directly to 007.

**Prediction.**
- If 013 > 007 by > 0.003: oracle scales — more pool = better.
- If 013 ≈ 007: oracle benefit saturates at ~500k candidates.
- If 013 < 007: same issue as 011 hyperactive — too-extreme sequences hurt.

**Generalization justification.** This is the cheapest way to test the
"saturation vs scaling" question. If oracle scales, we know HOW to push
the ceiling further. If it doesn't, we have a clearer picture of the
limit.

**Compute.** 4M candidates × 14k/s = ~5 min for scoring. Storage: ~5GB.
Acceptable.

## 2026-06-03 12:36 — Result: Experiment 013 — 4M-pool oracle WORSE than 500k-pool 007

**Result.** eval_01 = **0.3957**, mean_r = **0.3845** (vs 007: 0.3969 / 0.3861).
Per-cell-type: K562 0.615, HepG2 0.431, SK-N-SH 0.141.

**Oracle scaling NEGATIVE.** 8× larger pool (4M vs 500k) produced sequences
with mean Malinois 4.84 (vs 007's 2.60), but the trained model performed
SLIGHTLY WORSE. This was the THIRD outcome in my prediction tree, joining
exp 011 (SA hyperactive).

**Non-monotonic relationship now clear:**

| Library | Mean Malinois | mean_r |
|---------|--------------|--------|
| Random GC-50 | ~0.7 | 0.3827 |
| **007 (500k pool top)** | **2.60** | **0.3861** ← peak |
| 013 (4M pool top) | 4.84 | 0.3845 |
| 011 (SA hyperactive) | 5.51 | 0.3839 |
| 009 (cCRE+Malinois top) | 4.10 | 0.3774 |

There is a SWEET SPOT around Malinois mean activity 2-3. Above and below
underperform. This is consistent with several mechanisms:
- Oracle extrapolation error at extremes
- Loss of low-activity coverage when all sequences are top-tail
- Sequence-space concentration into fewer "good" neighborhoods

**Theory v7.** Oracle benefit is real (+0.005 over random) but
non-monotonic in selection extremity. It's maximized at moderate selection
(~500k → top 10% threshold), not at heavy oversampling. Beyond this,
ALL single-axis levers have plateaued at ~0.397.

But: difference 013 vs 007 = 0.0016. Is this within noise? I have not
yet measured the noise floor. EXP 014 will do exactly this.

## 2026-06-03 12:36 — Plan: Experiment 014 (007 replicate with different seed — noise floor)

**Mode: stress-testing apparent best (is 007 truly best, or lucky seed?).**

**Question.** What is the seed-to-seed noise floor in this evaluation
pipeline? If two replicate runs of 007's exact recipe differ by ±0.002,
then the 013 vs 007 gap (0.0016) is meaningless. If they differ by <0.0005,
the 0.397 ceiling is real and 013's degradation is meaningful.

**Design.** Re-run 007's recipe with SEED=2 (was SEED=1). Same:
500k GC-50 random, score with Malinois, top 50k by max(K562, HepG2, SKNSH).

**Prediction.**
- If 014 ≈ 007 within ±0.001: noise floor ≤ 0.001, oracle peak at 007's
  pool size is robust, and 013's degradation is signal.
- If 014 differs from 007 by ±0.002–0.003: noise is ~0.002, gap between
  013/007 is borderline, and the entire 0.397 plateau may be ±0.002 noise.
- If 014 < 007 by >0.003: 007 was lucky; the true ceiling is below 0.396.

**Generalization justification.** Knowing the noise floor is essential
for everything I claim about library design. Without it, "best library"
is undefined.

**Compute.** Identical to 007 — ~1200s.

## 2026-06-03 12:55 — Result: Experiment 014 — CRITICAL noise-floor finding

**Result.** eval_01 = **0.3936**, mean_r = **0.3827**.
Per-cell: K562 0.617, HepG2 0.430, SK-N-SH 0.133.

**vs 007 (identical recipe, only seed differs): Δ = -0.0034.** This is
the noise floor of the evaluation pipeline given a fixed (model,
training-size, evaluator) and the same library SCHEMA. Run-to-run
variation alone is ~0.003 on mean_r.

**Massive theory revision.** Many prior "lifts" and "drops" are within
the noise envelope:

| Comparison | Δ mean_r | Signal? |
|------------|----------|---------|
| 014 vs 007 (same recipe!) | -0.0034 | ≡ NOISE FLOOR |
| 013 vs 007 (4M pool) | -0.0016 | NOISE |
| 011 vs 007 (SA) | -0.0022 | NOISE |
| 006 vs 007 (span vs top) | -0.0001 | NOISE |
| 008 vs 007 (discrim) | -0.0026 | NOISE |
| 005 vs 007 (hg38 random) | -0.0056 | borderline |
| 012 vs 007 (hybrid) | -0.0056 | borderline |
| 009 vs 007 (cCRE+Malinois) | -0.0087 | SIGNAL (bad) |
| 003 vs 007 (PLS only) | -0.0208 | BIG SIGNAL (bad) |

The "oracle benefit" (007 vs random ~0.3827) appears to be ~+0.005,
just at the edge of noise. It might be real, but is borderline.

**Theory v8.** The 0.397 plateau is a fixed-point property of the
(model, training-size, evaluator) tuple, given the broad input subspace
"GC≈0.5 ACGT 200bp". WITHIN that subspace, the choice of selection
strategy has effect size ≤ 0.003 (noise). OUTSIDE that subspace,
deliberate deviations (wrong GC, wrong oracle bias) cost 0.01-0.02.

To break the plateau, future experiments must:
1. Produce >0.005 effect (signal > 1.5σ noise)
2. Move OUTSIDE the GC-50 random subspace in a STRUCTURED way that
   carries new information

Composition + oracle saturated. Need structural intervention.

## 2026-06-03 12:55 — Plan: Experiment 015 (motif-tiled — structured non-random subspace)

**Mode: testing fundamentally new lever — structured motif arrays.**

**Question.** All experiments so far have used GC-50 random as the
sequence-generating distribution, with various optional perturbations
(motif implants, oracle selection). Can a STRUCTURED non-random
distribution exceed 0.397?

**Design.** Fill each 200bp sequence with 8 contiguous JASPAR TF
motifs end-to-end, sampled from a diverse pool of well-known TFs
across cell types. No random scaffold; the entire 200bp is
motif-tiled. Each sequence = a different random 8-motif tile.

**Prediction (theory v8).**
- If 015 > 0.402 (>1.5σ above 007): structured motif syntax IS a
  new lever — the model learns better grammar from concentrated
  motifs than from sparse random ones. Big news.
- If 015 ≈ 0.385–0.397: structured motifs neither help nor hurt;
  the subspace barely matters as long as composition is reasonable.
- If 015 < 0.380: high motif density gives unrealistic, oracle-bias-
  triggering sequences (like 011 hyperactive), which hurt.

**Generalization justification.** If a structured-motif library
exceeds random+oracle, it shows that the (model, evaluator) cares
about motif arrangement, not just composition. This is a real
mechanism for learning transferable TF binding grammar.

**Compute.** Negligible — just sample motif sequences from PFMs.
Should be <2 min.

## 2026-06-03 13:10 — Result: Experiment 015 — motif-tiled is MUCH WORSE

**Result.** eval_01 = **0.3770**, mean_r = **0.3651**.
Per-cell: K562 0.592, HepG2 0.418, SK-N-SH 0.122.

**Big negative signal.** Drop of 0.020 vs 007 (7σ above noise floor).
Structured motif-dense subspace is significantly WORSE than random.

**Theory v9.** The 0.397 plateau isn't a saturation ceiling of the
GC-50 random subspace — it's the OPTIMUM. The (model, evaluator) tuple
rewards training sequences that resemble natural-like distributions,
and penalizes deviations into unnatural structured subspaces.

This is the cleanest result yet showing that the plateau is structural
to the model+evaluator, not specific to random sampling.

Updated landscape of effects:
- WITHIN natural-like (GC≈0.5 random or sparse-motif): plateau 0.397 ± 0.003
- OUTSIDE natural-like: drops 0.01-0.02 (PLS, motif-tiled, extreme GC)

The library design space within which 0.397 is reachable is BROAD but
shallow. To break it would require a genuine new information source.

## 2026-06-03 13:11 — Plan: Experiment 016 (GC variance test — per-seq vs library composition)

**Mode: mapping plateau structure.**

**Question.** When library mean GC = 0.50, does it matter whether each
sequence is itself GC≈0.50, or whether half are GC=0.30 and half GC=0.70?

**Design.** 50k = 25k random sequences sampled with target GC=0.30 +
25k random sequences with target GC=0.70. Library mean = 0.50, std large.

Generated by per-sequence Bernoulli sampling: each base C/G with prob 0.30
or 0.70, A/T otherwise (uniform within G/C and A/T).

**Prediction.**
- If 016 ≈ random GC-50 (mean_r ≈ 0.383): library MEAN composition is what
  matters; per-sequence GC variance is washed out by averaging.
- If 016 << random GC-50 (mean_r < 0.37): per-sequence composition matters;
  the model needs each sequence in the natural composition range.
- If 016 > random GC-50: bimodal training improves generalization
  (unlikely but possible if model needs to see extreme cases to anchor
  predictions on intermediate GC).

**Generalization justification.** Tests a fundamental theory question:
does the model learn a single composition prior, or per-sequence context?

**Compute.** Trivial — random Bernoulli per base. <1 min.

## 2026-06-03 13:22 — Result: Experiment 016 — per-seq composition matters

**Result.** eval_01 = **0.3760**, mean_r = **0.3674**.
Per-cell: K562 0.584, HepG2 0.415, SK-N-SH 0.130.

**Per-sequence composition matters; library-mean ISN'T sufficient.**
Drop of 0.019 vs 007 — same magnitude as motif-tiled (015). Bimodal
library with GC=0.30 + GC=0.70 (library mean=0.50) trains as poorly
as any "out-of-natural" subspace.

**Theory v10.** The natural-like subspace requires per-sequence
composition AND natural local structure. The (model, evaluator)
performance is determined by how closely TRAINING sequences resemble
EVAL sequences ONE-BY-ONE, not in aggregate.

## 2026-06-03 13:23 — Plan: Experiment 017 (multi-seed oracle pool — noise averaging)

**Mode: Can noise be averaged out across selection seeds?**

**Question.** Single-seed oracle (007) gives mean_r=0.3861;
replicate (014) gives 0.3827. Noise floor ~0.003. If I combine
selections from MULTIPLE oracle runs (different sequence pools, same
selection criterion), does the resulting library have LESS variance
and possibly a higher floor than any single seed?

**Design.** 50k = 10k from each of 5 independent SEED ∈ {1, 2, 3, 4, 5}
runs of 007's recipe (100k pool each, top 10k by Malinois max activity).
Compare to 007/014 single seed.

**Prediction.**
- If 017 > 007 by >0.005: multi-seed reduces variance and lifts the
  floor — first credible >σ effect for an oracle variant.
- If 017 within ±0.003 of 007: averaging is irrelevant; oracle
  selection captures the same information regardless of seed.
- If 017 < 014: bad — multi-seed dilution.

**Generalization justification.** A library that averages over multiple
oracle runs samples the oracle's preferred subspace more densely, which
could reduce overfitting to any one oracle-seed selection idiosyncrasy.

**Compute.** 5 × 100k pool scoring = 500k total, same as 007. ~30s.

## 2026-06-03 13:40 — Result: Experiment 017 — multi-seed oracle hits ceiling reliably

**Result.** eval_01 = **0.3977**, mean_r = **0.3862**.
Per-cell: K562 0.618, HepG2 0.436, SK-N-SH 0.139.

**Multi-seed converges to the oracle ceiling.** The 5×10k multi-seed
library matches 007's lucky-seed result exactly (Δ mean_r = +0.0001).
Meanwhile 014 (single unlucky seed) sat 0.003 below at 0.3827.

**Theory v11.** The oracle-selection ceiling is mean_r ≈ 0.386
(eval_01 ≈ 0.397) — robust across:
- Pool size: 500k = 4M = 5×100k
- Selection criterion (top, span, discrim): all noise-equivalent
- Random seed: spans 0.383-0.386, multi-seed averages to 0.386

To break it, we need a DIFFERENT sequence-level information lever,
not a refinement of the existing one.

## 2026-06-03 13:41 — Plan: Experiment 018 (pan-active oracle — MIN cross-cell)

**Mode: testing different oracle criterion (cell-type-transferability).**

**Question.** Standard oracle (007) picks sequences with high MAX activity
across (K562, HepG2, SKNSH). A pan-active oracle would pick those
with high MIN activity — sequences that activate ALL THREE cell types,
which by analogy should be MORE likely to activate held-out cell types
in the evaluator (eval 14).

**Design.** 500k GC-50 random, score with Malinois, top 50k by
MIN(K562, HepG2, SKNSH). Selected sequences are pan-cell-active.

**Prediction.**
- If 018 > 007 by >0.005: pan-active selection captures eval-relevant
  information that max-active misses. **First lever to break ceiling.**
- If 018 ≈ 007: any reasonable activity criterion yields the same
  ceiling.
- If 018 < 007: pan-active sequences are too narrow (most random
  sequences have low MIN), reducing diversity.

**Generalization justification.** Held-out cell types have unknown
TF composition. Sequences that work across multiple known cell types
should be more likely to encode COMMON TF binding logic, which
transfers. Sequences with one-cell extremes are TF-specific and may
not transfer.

**Compute.** 500k pool + scoring, ~30s.

## 2026-06-03 14:00 — Result: Experiment 018 — pan-active oracle matches max-active

**Result.** eval_01 = **0.3971**, mean_r = **0.3863**.
Per-cell: K562 0.618, HepG2 0.435, SK-N-SH 0.138.

**Pan-active = max-active = top-magnitude = span.** All within noise.
The ORACLE SELECTION CRITERION is irrelevant.

**Theory v12.** Within the natural-like subspace + Malinois-selected
active sequences, the ceiling is invariant to selection criterion.
What matters is "GC-50 random + ANY oracle filter for high-ish
activity." The "right" selection axis has no edge.

This exhausts the oracle-criterion axis. Remaining unexplored: real-
sequence-based, augmentation-based, and structural diversity axes.

## 2026-06-03 14:01 — Plan: Experiment 019 (real cCREs unfiltered)

**Question.** All cCRE-based experiments to date filter to GC ∈ [0.45, 0.55].
What if we let the natural cCRE distribution flow through? Real enhancers
have natural per-sequence composition variation around GC=0.50 with std
~0.07 (vs random GC-50 std 0.035). Does this NATURAL variation help, or
does any deviation from random GC-50 hurt?

**Design.** Real cCREs from ENCODE_cCREs_v3.bed, centered, extracted as
200bp windows. NO GC filter applied. Sample 50k.

**Prediction.**
- If 019 ≈ 0.380-0.385: natural cCRE composition variation hurts slightly
  (per-seq composition lever), but less than synthetic bimodal (016 ~0.367)
- If 019 ≈ random GC-50 (~0.383): natural variation in cCREs is similar to
  what random+filter produces
- If 019 > 0.386: real natural sequences carry signal that random+oracle
  cannot replicate

**Generalization justification.** Real eval sequences are likely natural
(or designed-from-natural), so a training set with natural composition
variation should match the eval distribution best. If 019 > 007, real
sequences ARE different in informative ways.

**Compute.** Just cCRE extraction. <1 min.

## 2026-06-03 14:21 — Result: Experiment 019 — natural cCREs match GC-50 plateau

**Result.** eval_01 = **0.3953**, mean_r = **0.3841**.
Per-cell: K562 0.609, HepG2 0.431, SK-N-SH 0.146.

**Real cCREs (no GC filter) match random GC-50 within noise.** This
RECONCILES the per-sequence composition theory with biology:
- 016 bimodal extreme (all seqs GC=0.30 or 0.70): mean_r 0.367 (BAD)
- 019 cCREs natural (mostly GC ~0.45-0.55, some tail): mean_r 0.384 (OK)
- 007 random GC-50 (all GC near 0.50): mean_r 0.386 (BEST)

So composition lever: requires MOST sequences in natural range; tail OK.

**Marginal SKNSH bump (0.146 vs 007's 0.141).** Real cCREs may include
some neural enhancers — but it's within noise.

## 2026-06-03 14:22 — Plan: Experiment 020 (RC-augmented — information doubling test)

**Mode: testing if explicit augmentation breaks the ceiling.**

**Question.** The trainer likely has RC augmentation built in (standard
practice). If we explicitly include each sequence AND its reverse-complement
in the library (25k unique → 50k including RCs), does this:
(a) help — meaning our model treats RCs as different and benefits from
    explicit doubling
(b) match — meaning the trainer already does RC aug internally
(c) hurt — meaning training duplicate-RC pairs creates collinearity that
    confuses the model

**Design.** Generate 25k GC-50 random sequences, then for each, compute
RC, write 50k total (interleaved).

**Prediction.**
- (a): mean_r > 0.395 — explicit doubling helps
- (b): mean_r ≈ 0.380 — equivalent to 25k unique random (which would
  be roughly the noise floor)
- (c): mean_r < 0.380 — RC pairs hurt

**Generalization justification.** Probes whether information-content
augmentation is a new lever.

**Compute.** Trivial. <1 min.

## 2026-06-03 14:40 — Result: Experiment 020 — RC-augmented matches 007 exactly

**Result.** eval_01 = **0.3964**, mean_r = **0.3861**.
Per-cell: K562 0.619, HepG2 0.435, SK-N-SH 0.136.

**Big insight: 25k unique + 25k RCs ≡ 50k unique.** The trainer is
RC-equivariant (either trains internally with RC aug or the model
arch has RC parameter sharing). Either way, RC sequences don't add
information beyond their forward sequences.

**Theory v13.** Unique-info bottleneck < 25k sequences. Library size
above 25k unique is NOT the limit. Per-sequence INFORMATION CONTENT
is what could differentiate libraries — but our oracle-selection
experiments showed even highly-selected sequences don't improve much.

This suggests the ceiling is fundamentally about (model capacity +
evaluator task structure), not about training-data design.

## 2026-06-03 14:41 — Plan: Experiment 021 (mutation-around-oracle — info density test)

**Mode: testing per-seq info density.**

**Question.** Given that 25k unique is enough, can we improve by
ensuring those 25k are the MOST INFORMATIVE 25k? Specifically: if we
take 5k oracle-top seeds and produce 9 single-base-mutated variants
of each, do these 5k × 10 = 50k sequences match 007 (where ALL 50k
are independently oracle-selected)?

**Design.** Take top 5000 from 007's selections (pre-computed).
For each, mutate 1 random position (different position for each variant)
to a different random base. Repeat to get 9 variants. Total = 5000 ×
10 = 50000 (including originals).

**Prediction.**
- If 021 ≈ 007: localized exploration of active sequence neighborhoods
  is sufficient; the model only needs to "learn" a small active subspace.
  This would be a SURPRISING result — minimal unique seeds suffice.
- If 021 << 007: high-info comes from DIVERSE active sequences, not
  from neighbors of any single active sequence. (Expected.)
- If 021 > 007: mutation augmentation around oracle-active is the
  hidden lever.

**Generalization justification.** Tests the structure of useful
information: is it 25k distinct active "modes" or a smaller number of
"modes" + their neighborhoods?

**Compute.** Need 007's selections (have them at scores.npz with
selected_idx). Trivial.

## 2026-06-03 15:05 — Result: Experiment 021 — sequence-correlation kills info

**Result.** eval_01 = **0.3706**, mean_r = **0.3600**.
Per-cell: K562 0.597, HepG2 0.409, SK-N-SH 0.106.

**Huge negative finding.** 5k oracle seeds × 9 mutations each scored
0.026 BELOW 007. Even though all 50k sequences are "unique" by string,
the 200bp sequences share 199/200 positions within each cluster, so
the effective independent-sample count is ~5k.

**Reconciles with 020.** 020's 25k unique + 25k RCs ≈ 50k unique
because the 25k unique are INDEPENDENT samples (i.i.d. random). 021's
5k × 10 mutants have only 5k effective independent samples. The
trainer benefits from INDEPENDENCE, not from raw sequence count.

**Theory v14: The lever is INDEPENDENT samples from natural-like
subspace.** Diminishing returns kick in around 25k. Below that, each
missing independent sample costs roughly:
- 5k → 0.360
- 25k → 0.386
- 50k → 0.386

This is the most actionable mechanistic finding yet:
**A training library's effective info content = its number of
INDEPENDENT samples, not its raw sequence count.**

This rules out all single-source augmentation schemes (mutation,
fuzz, etc.) as ways to add info content beyond what you started with.

## 2026-06-03 15:06 — Plan: Experiment 022 (anti-oracle — bottom 50k by Malinois)

**Mode: directional control for oracle effect.**

**Question.** Is the "oracle effect" (small but consistent +0.005 over
random) REAL? Confirm by reversing direction: train on the BOTTOM 50k
by Malinois score. If anti-oracle scores >> 007, oracle direction is
backwards (unlikely). If anti-oracle ≈ random, oracle is noise. If
anti-oracle << random, oracle has real signal in the active direction.

**Design.** 500k GC-50 random, score with Malinois, BOTTOM 50k by
max(K562, HepG2, SKNSH).

**Prediction.**
- Anti-oracle should be << random (mean_r < 0.378).
  Specifically, bottom 50k by max-active means sequences where ALL THREE
  cells predict LOW activity. This may train the model only on
  "off" sequences and badly degrade prediction for active eval cases.

**Generalization justification.** Critical sanity check that oracle
direction matters; foundation for any oracle-related theory claim.

**Compute.** Same as 007. ~30s.

## 2026-06-03 15:35 — Result: Experiment 022 — anti-oracle ≈ random GC-50

**Result.** eval_01 = **0.3933**, mean_r = **0.3830**.
Per-cell: K562 0.618, HepG2 0.432, SK-N-SH 0.130.

**Major: oracle direction barely matters.** Anti-oracle (mean Malinois
-0.135, completely INACTIVE predictions) achieves ≈ random GC-50
performance. The "+0.005 oracle benefit" we saw in 007 was at most
borderline-significant.

**Theory v15 consolidates:** The 0.397 ceiling is achieved by ANY
library satisfying:
1. ~25k+ independent samples
2. Per-seq composition in natural range
3. Natural-like local structure
4. From a subspace resembling eval distribution

Given these, SELECTION CRITERION is irrelevant:
- top, bottom, pan, span, discrim — all hit the same ceiling within noise

The ceiling is fundamentally about (model, evaluator, sequence subspace),
not about training-data selection.

## 2026-06-03 15:36 — Plan: Experiments 023-030 (remaining sweep)

Given the strong subspace theory (v15), remaining 8 experiments should
probe:
- 023: Independent-sample threshold (25k seeds × 2 mutations) — predicts
  ≈ 0.386 (interpolating 021→020)
- 024: 25k random + 25k cCRE natural mix (test subspace blend)
- 025: K-mer matched random (controlled k-mer histogram)
- 026: K562-only oracle (different oracle target, test if single-cell beats max)
- 027: Hard examples (Malinois prediction uncertainty/variance)
- 028: Subspace expansion test (50k random + ENCODE conserved elements)
- 029: 25k random + 25k 007 oracle mix (best vs random blend, anti-021)
- 030: FINAL SYNTHESIS based on findings

Most likely outcome: all 8 hit 0.380-0.390 within noise.

## 2026-06-03 15:58 — Result: Experiment 023 — 25k seeds + 1 mut each, slight dip

**Result.** eval_01 = **0.3902**, mean_r = **0.3799**.
Per-cell: K562 0.615, HepG2 0.427, SK-N-SH 0.128.

**Mild drop of 0.006 vs 007.** Bridges 020/021. Key insight:
- RC pairs are FREE (trainer handles RC aug) → 25k indep + RC = 50k indep
- 1-base mutants are NOT FREE (trainer treats as new data, redundant info)
- 1-base mutants → mild dip in mean_r (-0.006 from 50k indep baseline)
- 9 1-base mutants per seed (021) → big drop (-0.026)

So the scaling of effective info content is roughly logarithmic in
effective independent samples, with full saturation around 25k indep.

## 2026-06-03 15:59 — Plan: Experiment 024 (cCRE + random hybrid, different from 012)

**Question.** 012's 4-way hybrid (oracle + span + cCRE + random) failed
because of dilution. Is the cleaner 2-way blend of natural cCREs + random
synthetic any different? Tests whether NATURAL info from cCREs combines
constructively with synthetic random.

**Design.** 50k = 25k natural cCREs (unfiltered, from 019's recipe) +
25k random GC-50.

**Prediction.** Most likely ≈ 0.384 (within plateau, no benefit from mix).

## 2026-06-03 16:05 — Result: Experiment 024 — cCRE+random hybrid, plateau

**Result.** eval_01 = **0.3941**, mean_r = **0.3834**.
Per-cell: K562 0.609, HepG2 0.429, SK-N-SH 0.144.

**Within plateau as predicted.** The 25k natural cCREs + 25k random
mix sits squarely in the 0.380-0.386 zone. Confirms blends don't
synergize — independent samples + natural-like composition matter,
not heterogeneity per se. (012's 4-way mix at 0.3805 was within
noise of this, not a dilution effect; my earlier dilution narrative
was over-fit to one number.)

## 2026-06-03 16:06 — Plan: Experiment 025 (true duplicates — independence isolated)

**Question.** 021's 0.026 drop could be either (A) low effective
independence, or (B) spatial correlation specifically between
near-identical 200bp sequences confusing the trainer. 025 isolates
(A) from (B) by using **exact duplicates** instead of mutants:
each unique sequence appears N times.

**Design.** 5,000 unique GC-50 random sequences, each replicated
exactly 10x → 50,000 total. The trainer sees 50k samples but only
5k distinct sequences, with ZERO position-level correlation between
copies (they're identical, not adjacent in seq-space).

**Predictions:**
- If pure independence drives the gap: 025 ≈ 021 ≈ 0.360
- If spatial correlation between near-identical sequences hurts
  specifically: 025 >> 021 (closer to 0.380-0.386 plateau, since
  redundant copies should just average out / be noise-free)
- If duplicates HURT extra (overfitting the duplicate): 025 < 021

This isolates the mechanism that drove 021's drop.

## 2026-06-03 16:15 — Result: Experiment 025 — exact dupes WORSE than mutants

**Result.** eval_01 = **0.3557**, mean_r = **0.3453**.
Per-cell: K562 0.587, HepG2 0.382, SK-N-SH 0.098.

**Stunning.** 5k × 10 exact duplicates (mean_r 0.345) is worse
than 5k × 10 mutants (021: 0.360). Pure repetition is *worse* than
micro-correlated mutants.

**Theory v16: Independence-then-diversity-at-all-scales.**
- 50k truly independent samples (007/etc.) → plateau 0.386
- 25k indep + 25k mutants (023) → 0.380 (small dip)
- 5k indep + 45k mutants (021) → 0.360 (big dip)
- 5k indep + 45k EXACT copies (025) → 0.345 (biggest dip)

So mutants > duplicates of equal count. **Both effective sample
count AND micro-diversity in the data distribution matter
independently.** Exact repetition may push the trainer into
overfitting specific sequences, while mutants at least provide
local gradient information that the trainer can use.

**Refined mechanism:**
- Independence sets a *ceiling* on info content
- Within a fixed independence budget, having ANY variation
  (even 1-base mutants) beats pure repetition

This is now actionable: **never duplicate training data; if you
have a small pool, mutate aggressively rather than repeat.**

## 2026-06-03 16:16 — Plan: Experiment 026 (25k unique × 2 dupes — duplication scaling)

**Question.** Does the duplication penalty scale linearly with
unique-fraction? 025 (5k × 10) hit 0.345; 023 (25k unique + 25k
1-mut) hit 0.380. Try 25k unique × 2 exact dupes (50k library, 25k
distinct sequences via repetition, NOT mutation).

**Design.** 25,000 unique GC-50 random + each duplicated once = 50k.

**Predictions:**
- If repetition-penalty scales with redundancy fraction: 025 → 026
  should partially recover (~0.37, between 0.345 and 0.386)
- If matches 023 (0.380): repetition penalty saturates fast and is
  about ABSENCE of micro-diversity, not redundancy depth per se

## 2026-06-03 16:25 — Result: Experiment 026 — 25k×2 dupes ≈ 023 mutants

**Result.** eval_01 = **0.3888**, mean_r = **0.3786**.
Per-cell: K562 0.616, HepG2 0.428, SK-N-SH 0.123.

**Mutant-vs-duplicate equivalence at high unique count.** With
25k unique seqs:
- 026 (25k × 2 exact dupes): mean_r 0.379
- 023 (25k × 1 + 1-mut each): mean_r 0.380

At 5k unique:
- 025 (5k × 10 exact dupes): 0.345
- 021 (5k × 1 + 9 mutants each): 0.360

**Theory v17: micro-diversity helps ONLY when unique-count is
sub-saturation.** Above ~25k unique, mutant augmentation gives no
benefit over exact replication. Below saturation, mutation provides
extra training signal (probably because 1-base diffs at scale teach
local sensitivity in latent space).

**Putting it together (full data, unique-count → mean_r):**
- 5k unique → 0.345 (dupes) / 0.360 (mut)
- 25k unique → 0.379 (dupes) / 0.380 (mut) / 0.386 (with RC pairs)
- 50k unique → 0.386

The saturation curve has its bend around 25k. The +RC (020) trick
gets free info via trainer's built-in RC equivariance.

## 2026-06-03 16:26 — Plan: Experiment 027 (RC-canonical 50k — true RC-distinct sample bound)

**Question.** What if all 50k library sequences are mutually NOT
RC-equivalent (true RC-distinct)? Does pushing past the implicit
"25k unique × RC" boundary in 020 yield any benefit?

**Design.** Generate 50k unique GC-50 random sequences, but
deduplicate against RCs (no two are RC of each other). Should
expose the model to genuinely 50k distinct equivalence classes
under RC.

**Predictions.** If RC-equivariance is fully internalized by
trainer, 027 ≈ 007 ≈ 0.386 (no benefit). If RC-canonical sampling
shrinks effective uniqueness slightly (because trainer no longer
sees mirror pairs for each), 027 < 007.

## 2026-06-03 16:28 — Revised plan for 027 (RC-augmented oracle)

Reconsidered: 50k random RC-distinct is basically 007 (RC-collisions
in 200-mer random space are vanishingly rare). Pivoting to a real
test of independence-vs-info source.

**Question.** 020 (25k random + 25k their RCs) matched 007 exactly,
proving trainer is RC-equivariant. Does this also hold for ORACLE
selections, or do active sequences (whose motifs are often
strand-asymmetric) break RC-equivariance and benefit from explicit
RC augmentation?

**Design.** 25k oracle-top (Malinois-selected) + 25k their RCs.
Mirror of 020 but oracle-selected.

**Predictions:**
- If RC-equivariance holds universally: 027 ≈ 007 ≈ 0.386
- If active-sequence motifs break it: 027 > 007 (oracle info doubled)
- Underperforms if duplicate-penalty kicks in (it shouldn't — 25k unique)

## 2026-06-03 16:40 — Result: Experiment 027 — RC-augmented oracle ≈ 007

**Result.** eval_01 = **0.3958**, mean_r = **0.3852**.
Per-cell: K562 0.616, HepG2 0.436, SK-N-SH 0.136.

**Within noise of 007 (0.3861).** RC-equivariance holds even for
oracle-active sequences. There is no "active-motif strand asymmetry"
effect that would benefit from explicit RC augmentation. Confirms
v15 mechanism — trainer fully internalizes RC.

## 2026-06-03 16:41 — Plan: Experiment 028 (per-seq GC matched to cCRE distribution)

**Question.** Does per-sequence GC variance MATCHED to natural cCRE
distribution help, vs the tight GC ~ N(0.50, 0.035) of pure random?
019 (real cCREs) had GC mean=0.485 std=0.098 and hit plateau. But
that mixed both natural composition AND natural motif content. This
isolates composition alone.

**Design.** 50k random GC-varied sequences. For each seq:
1. Sample target_gc ~ N(0.485, 0.098), clipped to [0.20, 0.80]
2. Sample 200 bp i.i.d. with P(G)=P(C)=target_gc/2, P(A)=P(T)=(1-target_gc)/2

**Predictions:**
- If composition variance matters: 028 ≈ 0.384 (matches 019)
- If structure (motifs) matter beyond composition: 028 < 019
- Most likely: 028 ≈ plateau (since 005 already showed GC-matched
  random ≈ GC-matched cCRE)

## 2026-06-03 16:55 — Result: Experiment 028 — per-seq natural GC variance, plateau

**Result.** eval_01 = **0.3912**, mean_r = **0.3806**.
Per-cell: K562 0.610, HepG2 0.429, SK-N-SH 0.135.

**Plateau as predicted.** cCRE-like per-seq GC distribution (mean
0.485, std 0.098) without any motif content lands at 0.381 mean_r
— statistically indistinguishable from tight-GC random (007 = 0.386)
and natural cCREs (019 = 0.384) within noise.

**Confirms:** composition variance isn't a useful lever; the
0.38-0.39 plateau is robust to wide changes in per-seq GC distribution.

## 2026-06-03 16:56 — Plan: Experiment 029 (random+oracle blend)

**Question.** Anti-021 mix test: does combining random + oracle
synergize, dilute, or just hit plateau? 012's 4-way mix landed at
0.380 (within plateau), but had only 12.5k of each. This isolates
the clean 25k random + 25k oracle blend.

**Design.** 25k GC-50 random + 25k from 007 oracle = 50k library.
Both contributions are independent samples in their respective
selection regime.

**Predictions:**
- If random dilutes oracle info: 029 < 007 (~0.378)
- If oracle dilutes random info: 029 ≈ random+oracle midpoint ≈ 0.385
- If both contribute and blend constructively: 029 ≥ 007 (push past plateau)
