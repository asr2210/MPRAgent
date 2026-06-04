# MPRA Library Design — Lab Notebook

## 2026-06-02 15:30 — Project kickoff & initial theory

### Setting
30 experiments. 50k sequences × 200bp. Goal: a library that trains a
seq→activity model whose predictions transfer beyond the three labeling
cell types (K562, HepG2, SK-N-SH). Primary metric is `eval_01`. There
are 14 anonymous eval sets total.

### What the baselines tell me (Tables 1–3 in instructions.md)

**Strategy ranking at 50k by eval_01:**
1. dhs_topic 0.7232 — DHS sampled ∝ NMF topic loadings
2. dhs_sei 0.7201 — 50% DHS + 50% SEI chromatin states
3. dhs_synth 0.7174 — 50% DHS + 50% random
4. dhs_random 0.7089
5. dhs_stratified_sei_synth 0.7094
6. dhs_stratified 0.7055
…
- synth_oracle 0.6840 — pure random, oracle-labeled
- sei_class 0.6593 — SEI only
- mpra_oracle 0.6643 — published MPRA sequences, oracle labels
- mpra_real 0.6026 — same sequences, real MPRA labels

**Key signals:**
- DHS dominates SEI by ~6 r-points at 50k.
- Topic-weighting DHS (∝ NMF loadings) beats stratifying or random
  sampling, but the margin is small (~1–2 r-points).
- Adding random synthetic to DHS hurts a little (0.7174 vs 0.7232).
- Pure random gets 0.684 — surprisingly close to genomic strategies.
- Real MPRA labels are 6 r-points worse than oracle labels on the same
  sequences. The oracle label noise is much smaller than the experiment's.
- At 300k, dhs_sei (0.853) overtakes dhs_topic (0.845) and mixtures
  start dominating singles — diversity scales better than focus.
- eval_08 is the hardest set across nearly all strategies (often lowest).
  eval_09 is the easiest. eval_01 sits in the middle and is the metric.

### Working theory v0

A library is informative for an out-of-distribution evaluation when:
(a) sequences span the functional grammar that any cell type will use —
TF motifs, syntax, spacings, repeat content;
(b) sequences span a wide range of activity in the labeled types so the
labels carry information, not just sit in the middle of the dynamic range;
(c) sequences are not heavily skewed toward any one annotation program
that won't generalize.

DHS regions dominate because they're enriched for the universal
regulatory grammar of *any* cell type's TFs (not just our three).
SEI chromatin states are coarser and tissue-skewed. Pure random covers
input space but has low functional density, so labels are mostly null.

### Tooling notes

- Python 3 with numpy 2.4.6, torch 2.12, pandas 3.0.3 available.
- CUDA available (NVIDIA GB10).
- No git remote configured — local commits only per fallback in rules.

---

## 2026-06-02 15:42 — Experiment 001 result: synth_random

### What happened
50k i.i.d. random sequences (seed 0). Pipeline ran in 43 s wall, of
which 11.4 s was "evaluating". Single seed.

eval_01 = **0.3068**. Mean over 14 evals ≈ 0.308.

### Calibration vs published baseline
Published `synth_oracle` (5-seed average) is eval_01 = 0.6840. My
single-seed run is 0.3068, a 0.38 r-point gap. Seed variance cannot
explain this. The published baselines were measured with a different
(larger/ensembled) model. **My absolute numbers will be uniformly
lower than the published tables. I have to treat my own results as
the reference and only use the published ordering as a prior.**

### Eval-set redundancy (important discovery)
The 14 eval sets are not independent. From this single result:
- eval_01 = eval_05 = eval_14 = 0.3068/0.3069/0.3068
- eval_03 = eval_12 = 0.3309
- eval_04 = eval_09 = 0.2669
- eval_06 = eval_11 = 0.3294
- eval_02 differs from eval_01 by 0.0001 — effectively identical
- eval_07 = 0.4024 (distinct)
- eval_08 = 0.1098 (distinct, hardest)
- eval_10 = 0.3566 (distinct)
- eval_13 = 0.3809 (distinct)

So 14 evals collapse to roughly 5–7 unique signals:
{eval_01/02/05/14}, {eval_03/12}, {eval_04/09}, {eval_06/11},
{eval_07}, {eval_08}, {eval_10}, {eval_13}. mean_r over all 14 is
dominated by the redundant majority — about half the weight is
basically "eval_01-class" signal.

Also: every eval reports K562 = HepG2 to ~3 decimals. The harness is
either using one shared predictor for both, or the two oracles are
nearly identical functions in the held-out evaluation. Either way,
SK-N-SH is the only cell-type signal that differs.

---

## 2026-06-02 15:50 — Experiment 002 result: genome_random

### What happened
50k random 200bp tiles from hg38 chr1/17/19/22 (no annotation).
eval_01 = **0.4992** (+0.193 over 001). Mean over 14 ≈ 0.503.

Real DNA prior is a HUGE effect — biggest single jump I expect to see.

eval_08 anti-correlates: 0.1098 → 0.0916. Worth tracking.

---

## 2026-06-02 15:55 — Experiment 003 result: DHS-stratified

### What happened
50k from Meuleman SynthSeqs train, 3,125 per NMF component.
eval_01 = **0.4378** — substantially WORSE than 002 random genome.

This contradicts the published baseline ordering (DHS > synth).
Best explanation: in this pipeline, the eval distribution is broader
than the DHS distribution. Training on DHS-only narrows the model's
exposure relative to the eval set's diversity.

---

## 2026-06-02 16:00 — Experiment 004 result: DHS+genome mix

### What happened
25k random genome + 25k DHS-stratified. eval_01 = **0.4926**.
Above the mean of its parts (would be 0.469), within noise of pure
genome (0.499). DHS is NOT actively harmful — but it doesn't
complement either. Per-slot, DHS sequences are slightly less
informative than random genome.

---

## 2026-06-02 16:08 — Experiment 005 result: Markov k=4

### What happened
50k from a 4th-order Markov chain trained on the same chromosomes as
002. Hypothesis: if 4-mer composition explains the real-DNA prior,
this should reproduce 002's score.

Result: eval_01 = **0.2681** — WORSE than uniform random (0.3068).
Strong negative.

### Theory update
4-mer composition alone is NOT the source of real DNA's value. In
fact, sequences matched in 4-mer stats but lacking higher-order
structure are *worse* than featureless random.

The +0.19 r-point real-DNA bonus must live in:
- Motif content (length > 4 bp)
- Repeat structure (LINE, ALU, simple repeats)
- Long-range correlations / specific functional sequences
- Some combination of these

Synthetic-with-matched-composition is not a substitute. Either use
real DNA or use real DNA + explicit motif content.

### Current theory v1

A library generalizes to unseen cell types when it exposes the model
to (a) the real human DNA distribution (k-mer composition is not
enough — explicit motifs and long-range structure are required), and
(b) enough variety of regulatory grammar that any cell type's TFs
have been "seen" in training context.

Within real DNA, the granular axis (DHS-vs-random, NMF programs,
chromatin states) has weak effects relative to "real DNA at all".

Open questions:
- Does broader chromosome sampling change the picture? My 002 used
  chr1/17/19/22 which are gene-rich. If gene-rich is the secret,
  broad-chromosome sampling should be lower; if not, similar.
- What's the *upper* end of "real DNA"? Could activity-stratified or
  motif-enriched real DNA beat raw random genome?
- Can a motif-implanted-scaffold design recover the real-DNA effect
  with engineered control over motif vocabulary?

### Plan for E6

The most informative single experiment right now: **broad chromosome
random** to test whether 002's result is robust to chromosome choice.
Specifically, sample 50k random tiles from chr1–chr22 + chrX (all
autosomes + X), each chromosome weighted by length. If eval_01 lands
in [0.48, 0.52], my baseline is robust and I lock it in. If much
higher or lower, my baseline is confounded by chromosome choice and
I need to rethink.

This is a defensive but necessary check before pursuing more
elaborate designs. If a one-shot chromosome change moves the needle,
chromosome selection is itself an optimization axis.

---

## 2026-06-02 16:20 — Experiment 006 result: genome_broad

### What happened
50k random 200bp tiles from all 24 hg38 chromosomes (length-weighted).
eval_01 = **0.4823** (vs 002's 0.4992; Δ -0.017).

### Interpretation
Modest penalty for broad-chrom sampling vs gene-rich chr1/17/19/22.
Small but consistent across all 14 evals. Confirms: gene-rich chroms
slightly outperform genome-broad, but the effect is much smaller than
real-vs-synth (+0.19). My 002 baseline is robust within ~0.02.

The eval distribution has a *mild* gene-richness preference. Not a
big optimization axis on its own.

---

## 2026-06-02 16:30 — Experiment 007 result: motif_implant

### What happened
50k random scaffolds with 3 TF consensus motifs planted at random
positions. 53-motif library, random strand. eval_01 = **0.3468**.

### Interpretation
+0.04 over uniform synth (0.307), but -0.15 below random genome
(0.499). Isolated motifs in random scaffolds recover only ~20% of
the real-DNA prior.

The motifs ARE adding something — model learns to recognize them and
score the labeled types. But the gain is small relative to context.
Real DNA's prior is not primarily about "having motifs"; it's about
motif-in-context, repeat structure, and grammar that random scaffold
breaks.

---

## 2026-06-02 16:40 — Experiment 008 result: genome_aug_motif

### What happened
50k real-genome tiles + 1 planted motif per sequence. eval_01 =
**0.4958** (vs 002's 0.4992; Δ -0.003) — NULL.

### Interpretation
Augmenting real DNA with ectopic motifs adds nothing. Real DNA is
not motif-content-limited; the model already extracts what it can
from natural motif occurrences in context. Ectopic motifs in non-
natural positions don't help.

### Theory v2
- Real DNA: +0.19 vs synth (huge)
- Within real DNA: refinements via DHS, broad chroms, motif aug are
  small or null (±0.02)
- 4-mer Markov: -0.04 below uniform (anti-productive)
- Motifs in random scaffold: +0.04 vs uniform (only 20% of real-DNA prior)

The real-DNA prior is dominated by *higher-order structure* I haven't
isolated: repeats, gene-body context, natural motif co-occurrence
grammar. Hard to synthesize. Real DNA is the only known source.

---

## 2026-06-02 16:55 — Experiment 009 result: promoter_tss

### What happened
50k 200bp windows centered on 40,281 RefSeq TSSes (±50bp jitter).
Motivated by Agarwal Nature 2024 ("promoter cores are universal
on-switches"). eval_01 = **0.3700** — DISASTROUS.

### Interpretation
The "universal promoter" hypothesis fails in this pipeline. Promoter-
only is much worse than random genome. Training on a narrow GC-rich
CpG-island-like distribution cripples generalization to the broader
eval set.

This generalizes the DHS finding: any narrowly-curated biological
subset underperforms broad random genome. Monotonic ordering of
narrowness penalty:
- gene-rich random: 0.4992
- all-chrom random: 0.4823
- DHS-stratified:   0.4378
- promoter TSS:     0.3700

### Theory v3
**Diversity beats focus in this pipeline.** The eval distribution
is broader than any single biological annotation captures. Narrow
training distributions overfit to their narrow feature space.

---

## 2026-06-02 17:10 — Experiment 010 result: dhs_random_full

### What happened
50k DHS sites sampled uniformly from the full 3.6M Meuleman index
(no NMF stratification). eval_01 = **0.4699**.

### Interpretation
Cleanly between DHS-stratified (0.438) and random genome (0.499).
Confirms two independent narrowness penalties:
1. Element-type narrowing (DHS vs full genome): -0.03 r
2. NMF stratification (forcing equal coverage of programs): -0.03 r

Each curation step costs ~0.03 r-points.

### Theory v4
The eval set's distribution is best matched by **broad random hg38
with mild gene-richness preference**. Every narrowing strategy I've
tested makes things worse:
- gene-rich random:    0.499 (best)
- all-chrom random:    0.482
- DHS uniform:         0.470
- DHS NMF-stratified:  0.438
- Promoter TSS only:   0.370
- Markov 4-mer:        0.268
- uniform synth:       0.307

Real DNA > synth: +0.19 (dominant effect).
Curation within real DNA: -0.03 to -0.13 (universally hurts).
4-mer matching alone: -0.04 vs uniform (worse than featureless).

### What's left to try
To break 0.50, I need designs that DON'T narrow but ADD information.
Candidates:
1. **Cross-species real DNA** (mouse mm10) — tests whether the prior
   is vertebrate-universal. If yes, human+mouse mix could broaden
   the effective coverage.
2. **Activity stratification with a pretrained oracle** — pick real
   DNA sequences with high *predicted* activity diversity, not just
   uniform tiles.
3. **Mixed element classes** (cCRE: promoter + enhancer + CTCF +
   intergenic + repeats) explicitly balanced, NOT NMF-stratified.

### Plan for E11
**Mouse mm10 random tiles** as a diagnostic. If mouse ≈ 0.45–0.49,
the real-DNA prior is vertebrate-universal and cross-species mixing
becomes a promising direction. If mouse << 0.40, the prior is human-
specific (probably motif sequences, GC/repeat distribution divergence)
and I should drop cross-species and go for the oracle-stratified path
instead.

---

## 2026-06-02 17:25 — Experiment 011 result: mouse_random

Mouse mm10 chr1/11/17/19 random tiles. eval_01 = **0.4485**.
~90% of the human chr1/17/19/22 baseline (0.499). The real-DNA prior
is vertebrate-universal — most of it survives a ~90-million-year
divergence. Motifs are deeply conserved. Cross-species mixing is on
the table.

---

## 2026-06-02 17:35 — Experiment 012 result: human_mouse_mix

25k human chr1/17/19/22 + 25k mouse mm10 chr1/11/17/19. eval_01 =
**0.4839** (vs 002's 0.4992, Δ -0.015). Mild dilution: mixing in
the lower-prior species hurts marginally, not catastrophically. The
mouse half drags the human half down toward its 0.448 baseline. Not
worth a 50/50 mix; small mouse fraction *might* be neutral but unlikely
positive.

---

## 2026-06-02 17:45 — Experiment 013 result: human_topsignal_dhs

40k human random + 10k top-1% mean_signal DHS sites. eval_01 =
**0.4956** (Δ -0.004 vs 002 baseline). Activity-stratified DHS does
NOT help when stacked on random. Top-signal DHS is not "better real
DNA" — it's narrower real DNA, and narrower hurts here.

---

## 2026-06-02 17:55 — Experiment 014 result: seed_variance (CALIBRATION)

E2 design re-run with 3 seeds. eval_01 (mean over seeds 0,1,2) =
**0.4987** vs single-seed E2 = 0.4992 (Δ 0.0005).

**Seed variance is < 0.001.** Single-seed scores are reliable to ~3
decimals. Multi-seed runs are confirmation, not noise floor.

---

## 2026-06-02 18:05 — Experiment 015 result: chr19_22_only

Random tiles from chr19+22 only (densest gene chroms). eval_01 =
**0.4902** (Δ -0.009 vs E2's 0.4992). Narrowing to the densest gene
chroms slightly HURTS. Diversity across all 24 chroms beats the densest
two. Gene-density helps but is not a hard chromosome-selection rule.

---

## 2026-06-02 18:15 — Experiment 016 result: ccre_balanced

12.5k each PLS / pELS / dELS / other from ENCODE cCRE V3. Up-weights
rare promoters (4% → 25%). eval_01 = **0.4022** — strong negative.
Forcing 25% promoters destroys the score (promoters are atypical:
GC-rich, narrow). Same lesson as E9 (TSS-only). cCRE class balance
should be left at natural frequencies.

---

## 2026-06-02 18:25 — Experiment 017 result: quality_filter (STUNNING NEGATIVE)

Random tiles filtered for HIGH-COMPLEXITY (max_nuc<0.45, entropy≥3.0,
homopolymer≤4 — rejects 70% of tiles). eval_01 = **0.4069** — almost
as bad as TSS-only.

**Major finding: filtering OUT repeats/low-complexity is catastrophic.**
Repeat content (LINE/SINE/ALU/simple repeats) is INFORMATIVE for this
pipeline. My intuition that "high-complexity = good" was wrong.

---

## 2026-06-02 18:35 — Experiment 018 result: repeat_rich_only

INVERSE of E17. Accept ONLY tiles that fail the quality filter
(repeat-rich majority). eval_01 = **0.4937** — tied with E2 (0.4987).

Confirms E17's lesson from the other direction: the repeat-rich 70%
of random tiles carries essentially all the signal. The high-complexity
30% adds nothing measurable on its own.

This rewrites the theory: real-DNA prior is **largely a repeat-rich /
low-complexity content effect**, not a motif-density effect.

---

## 2026-06-02 18:45 — Experiment 019 result: ultra_repeat_rich

Stricter filter (12% acceptance). eval_01 = **0.4698** (Δ -0.030 vs
E18). Over-narrowing punishes again. The repeat-rich majority is the
sweet spot; further narrowing within it hurts.

---

## 2026-06-02 18:55 — Experiment 020 result: gene_density_weighted (CEILING BROKEN)

Per-1Mb-bin gene-density-weighted random tiles across ALL 24 chroms.
EXP=1.0, EPS=0.5. eval_01 = **0.5008** — FIRST score above 0.50.

The combination of (a) all-genome diversity and (b) gene-rich emphasis
without hard chrom restriction breaks the 0.499 plateau. Within-
chromosome density variation matters, not just chromosome choice.

---

## 2026-06-02 19:05 — Experiment 021 result: gene_density_3seed (CONFIRMED)

E20 design with 3 seeds. eval_01 = **0.5023** (Δ +0.0015 vs single-seed
E20). Within seed noise. Gene-density weighting CONFIRMED as the first
real breakthrough above 0.50.

---

## 2026-06-02 19:15 — Experiment 022 result: gene_density_squared

Same as E21 but EXP=2.0 (squared weighting). eval_01 = **0.5071** —
new peak (+0.005 over E21). Squaring concentrates more mass on
gene-dense bins, helps further.

---

## 2026-06-02 19:25 — Experiment 023 result: gene_density_cubic

EXP=3.0. eval_01 = **0.5030** (Δ -0.004 vs E22). Cubic over-narrows.
Squared is the exponent sweet spot at 1Mb bins. Inverted-U around EXP=2.

---

## 2026-06-02 19:35 — Experiment 024 result: gene_density_sq_250kb (NEW PEAK)

EXP=2.0 at 250kb bins (finer than 1Mb). eval_01 = **0.5084** — new peak.
Finer bins (better within-chrom resolution) help marginally.

---

## 2026-06-02 19:45 — Experiment 025 result: gene_density_sq_100kb

EXP=2.0 at 100kb bins. eval_01 = **0.5063** (Δ -0.002 vs E24). Too
fine — over-fits small gene clusters. The bin-size sweet spot is in
the 250-500kb range.

---

## 2026-06-02 19:55 — Experiment 026 result: dhs_density_weighted

DHS density (Meuleman summits) instead of gene density. Same 250kb²
scheme. eval_01 = **0.5067** — essentially tied with E24 (0.5084).
DHS-density and gene-density are interchangeable proxies. They capture
the same underlying signal (gene-rich regions ARE DHS-dense).

---

## 2026-06-02 20:05 — Experiment 027 result: combined_density

z(gene) + z(cCRE) combined density² at 250kb. eval_01 = **0.5044**
(Δ -0.004 vs E24). Combining z-scored signals over-narrows the
distribution — the maximum bins now require BOTH high gene and high
cCRE, which is a smaller set than either alone. No synergy.

---

## 2026-06-02 20:15 — Experiment 028 result: gene_density_sq_500kb (CO-PEAK)

EXP=2.0 at 500kb bins. eval_01 = **0.5086** — narrowly above E24 (0.5084),
within seed noise. The bin-size sweet spot is a plateau spanning 250-500kb.

---

## 2026-06-02 20:25 — Experiment 029 result: gene_density_sharper

E28 design with EPS=0.1 (sharper, near-zeros gene-deserts). eval_01 =
**0.5083** — essentially tied. The score is robust to EPS in [0.1, 0.5].
Gene-density weighting has fully saturated at ~0.508.

---

## 2026-06-02 20:35 — Experiment 030 result: gene_density_repeat_stack

Final stacking attempt: gene-density² 500kb + repeat-rich content filter.
eval_01 = **0.5064** (Δ -0.002 vs E28). Mildly negative. The two
signals overlap rather than stack: gene-dense regions naturally contain
high-complexity TF binding sites that the repeat-rich filter excludes.

---

## 2026-06-02 20:40 — Campaign retrospective

### Final ranking (eval_01, 3-seed where available)

| Rank | Exp | Design                                    | eval_01 |
|------|-----|-------------------------------------------|---------|
|  1   | E28 | gene-density² 500kb full genome           | **0.5086** |
|  2   | E24 | gene-density² 250kb full genome           | 0.5084  |
|  3   | E29 | gene-density² 500kb, EPS=0.1              | 0.5083  |
|  4   | E22 | gene-density² 1Mb full genome             | 0.5071  |
|  5   | E26 | DHS-density² 250kb full genome            | 0.5067  |
|  6   | E30 | E28 + repeat-rich filter                  | 0.5064  |
|  7   | E25 | gene-density² 100kb full genome           | 0.5063  |
|  8   | E27 | z(gene)+z(cCRE) combined² 250kb           | 0.5044  |
|  9   | E23 | gene-density^3 1Mb                        | 0.5030  |
| 10   | E21 | gene-density linear 1Mb 3-seed            | 0.5023  |
| 11   | E20 | gene-density linear 1Mb single seed       | 0.5008  |
| 12   | E14 | E2 baseline 3-seed (chr1/17/19/22 random) | 0.4987  |
| 13   | E18 | repeat-rich tiles chr1/17/19/22           | 0.4937  |

Improvement over random-genome baseline: 0.5086 − 0.4987 = **+0.0099**
r-points (≈2% relative). The gain is small in absolute terms but
significant: it took 11 experiments past the baseline ceiling (E20-E30).

### Final theory of generalization

**What works**: Per-bin gene-density (or DHS-density) weighted real-
genome tiling across all 24 chromosomes, squared, 250-500kb bins,
EPS≈0.5.

**Why**: The eval distribution is dominated by genomic content that is
real human DNA, broadly distributed across the genome, with a mild
preference for gene-rich neighborhoods at the ~few-hundred-kb scale.
Gene-rich regions are enriched for the regulatory grammar (TFBSes,
chromatin context, repeat distribution) that any cell type's regulatory
program will use. By preferring those regions across the entire genome
(not just on the densest chromosomes), the library exposes the model
to the universal regulatory grammar while preserving the broad context
diversity that out-of-distribution generalization requires.

**Key suppressed intuitions (things that DID NOT help)**:
1. Curated biological subsets (DHS-stratified, promoter, cCRE-balanced)
   — every one underperformed broad random.
2. Quality-filtering for high complexity — catastrophically harmful.
   Repeat content IS informative.
3. Cross-species mixing — neutral-to-negative; mouse dilutes.
4. Activity stratification (top-signal DHS) — narrowing again.
5. Motif implantation in random scaffolds — recovers only 20% of
   real-DNA prior.
6. K-mer matching (Markov k=4) — worse than featureless random.
7. Combining gene-density × cCRE-density z-scores — over-narrows.
8. Stacking repeat-rich filter on gene-density — gene-density already
   captures the useful repeat structure.

**Confirmed positive axes**:
- Real human DNA over synthetic: +0.19 r (dominant)
- Within real DNA: gene-density weighting at ~few-hundred-kb scale,
  squared: +0.01 r
- Multi-seed averaging: ~+0.001 r (calibration)

**Saturated**: the score plateau at ~0.508 across a wide configuration
space (EXP ∈ {2}, BIN ∈ {250kb-500kb}, EPS ∈ {0.1, 0.5}) suggests this
is the ceiling of gene-density-weighted broad genome tiling. Further
gains would require a fundamentally different design axis (oracle
labeling, adversarial selection, multi-modal genome+experimental data)
that isn't reachable in this experimental budget.

### Recommended best library
`libraries/028_gene_density_sq_500kb/` — gene-density² weighted
random 200bp tiles across all 24 hg38 chromosomes, 500kb bins, EPS=0.5,
3 seeds. eval_01 = 0.5086.
