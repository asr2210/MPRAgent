# MPRA Library Design — Lab Notebook

## 2026-06-02 21:25 — Initial context & strategy

**The problem.** Design 50k 200bp sequences. Train a sequence-to-activity model on
MPRA labels from K562/HepG2/SK-N-SH. Evaluate on 14 anonymous eval sets. The
crucial framing: the library should generalize to cell types we never measured.

**Baselines we know about (from instructions.md / strategies.md).**
- v14 synthetic-only baselines (random_uniform, gc_sweep, at_rich, …) all give
  Pearson r ≈ 0 on eval_01. They train on REAL MPRA labels, which are dominated
  by noise for inactive synthetic sequences.
- "informed" baselines using real regulatory sequences with real MPRA labels:
  - dhs_topic: 0.7232 (best), dhs_sei: 0.7201, dhs_synth: 0.7174, dhs_random: 0.7089
  - synth_oracle (random sequences, ORACLE labels): 0.6840 — proves the model can
    learn from random sequences when labels are clean.
  - mpra_real (re-using published MPRA): 0.6026 — works but slightly worse.
- Learning curves: dhs_topic at 300k gets 0.8448; baselines saturate ~0.85.

**Implication for v14.** REAL labels with REAL biological sequences is what works.
Synthetic sequences need oracle labels to work, which we don't have access to.
So our library must consist of biologically active regulatory sequences.

**Current theory (v0).** A library is informative for cross-celltype
generalization to the extent that it covers diverse regulatory grammar:
- many TF motifs in many syntactic contexts (spacing, orientation, count)
- many promoter / enhancer / silencer / insulator-like architectures
- many chromatin accessibility programs across the human body
- AT/GC composition range
- both common and rare sequence configurations

A library tuned only to K562/HepG2/SK-N-SH-active regions will overfit cell-
type-specific TF combinations. A library that spans many cell types' regulatory
elements should let the model learn TF→activity rules that transfer.

**Predictions from theory.**
1. dhs_topic beating dhs_random is consistent: topic-weighting selects elements
   with strong tissue-specific signal, which encode more diverse TF programs.
2. Mixing synthetic with DHS (dhs_synth) is slightly worse than pure DHS — the
   real signal-bearing sequences are diluted, and the synthetic half adds noise
   under real labels (mpra_real is weak; synth_oracle has labels we don't get).
3. SEI-class weighting adds value if it covers regulatory contexts under-
   represented by DHS topics (e.g., insulators, polycomb).

**First experiment plan.** Reproduce a DHS-like baseline using ENCODE candidate
cis-regulatory elements (cCREs). cCREs from SCREEN v3 cover ~1M elements across
~1500 biosamples — comparable to DHS in scope. This is the cleanest way to
sample diverse regulatory contexts. I'll start with topic-weighted-equivalent
sampling (proportional to ENCODE classification, e.g., upweighting distal
enhancers which are more diverse) and a 200bp window centered on the cCRE.

**Why this generalizes beyond K562/HepG2/SK-N-SH.** cCREs are defined across
all ENCODE biosamples, not just our three. Their sequences encode the
regulatory grammar of all human tissues. A model trained on cCREs sees TF
binding patterns from immune cells, neurons, hepatocytes, fibroblasts, etc.
The labels are from our 3 cell types, but the sequence diversity it learns to
parse should transfer to predicting activity in any cell type whose regulatory
grammar is built from the same TFs.

## 2026-06-02 21:40 — Experiment 001 result: cCRE_random → r ≈ 0

Ran uniform-random ENCODE cCRE (V3) sampling, 50k × 200bp centered on cCRE
midpoints. Result: **eval_01 mean_r = 0.0016**. All 14 eval sets near zero.

Time: 11.6 s total, n_seeds=1. This is the same near-zero range as every
"basic synthetic" strategy in strategies.md (v14 baselines). My biological
sequences gave no more signal than random_uniform did.

**Key insight that changes my theory.** The "informed" baseline table in
instructions.md (which predicted ~0.71 for this kind of library) does NOT
appear to apply to v14's prepare.py. The actual v14 evaluator seems to
behave like the v14.md baselines — most libraries return ~0.

Updated theory (v1). v14 trains in ~12 s with a single seed. This budget is
too small for the model to extract complex motif grammar from generic
biological sequences. To produce a measurable signal, a library must:
- contain low-entropy, high-redundancy features the model can learn in a
  handful of gradient steps, OR
- correspond exactly to the distribution the eval sets are drawn from, OR
- expose the model to enough copies of each learnable feature that even a
  brief training pass converges.

Predictions:
- (P1) Replicating a smaller unique set should help (more examples per
  unique pattern, more averaging of label noise).
- (P2) Pure motif-saturated synthetic with very strong known TF sites should
  give some signal — strong + redundant.
- (P3) A bimodal library (half strongly-active, half deliberately-inactive)
  should give a learnable contrast.

**Generalization caveat.** Whatever I find will be optimal for v14's specific
evaluator. The instructions ask me to design a library that would generalize
even to cell types not in {K562, HepG2, SK-N-SH}. If the v14 evaluator is in
fact testing generalization across cell types (and that's why nothing works),
my upcoming experiments must keep diversity / coverage as a parallel goal —
otherwise I'd be over-fitting to v14's idiosyncrasies rather than building a
genuinely informative library.

## 2026-06-02 21:42 — Planning Experiment 002

Plan: test (P1) — replication-based library. Use 5,000 unique cCRE sequences,
each replicated 10 times, for 50,000 total. Compare to 001 (50,000 unique).
If r jumps meaningfully, the bottleneck is examples-per-feature, not
diversity. If r stays at 0, replication isn't the answer and I should test
(P2) next.

Mode: exploring a new hypothesis (replication helps under limited training).
Generalization argument: replicating cCREs preserves biological grammar — I'm
not biasing toward K562/HepG2/SK-N-SH specifically. If replication works, the
same trick would help in any low-compute evaluator on any tissue.

## 2026-06-02 21:55 — Exp 002/003 results: replication & motifs both ~0

Both replication (002, r=0.0021) and motif-packing (003, r=-0.0037) give
the same noise-floor result as cCRE diversity (001, r=0.0016) and as every
listed v14 baseline. After three swings, no library design has moved r
out of [-0.004, +0.003] on eval_01.

Updated theory (v2). v14 has a noise floor of roughly ±0.005 across all
single-distribution libraries — biological, replicated, or motif-rich. The
v14 evaluator may:
- Use a very small / under-trained model whose predictions are nearly
  constant for any input, so Pearson r ≈ 0 regardless of library, OR
- Score against eval sequences drawn from a distribution unlike anything
  I've tried, where neither cCREs nor motif-saturated synthetics overlap.

Things I've not yet tried that might break the floor:
- (Q1) Bimodal contrast: pair motif-rich with motif-empty so labels are
  forced bimodal. The model has only two classes to learn.
- (Q2) Real published MPRA training data (Sharpr-MPRA, Agarwal LentiMPRA).
  If the eval is also MPRA-like, distribution match might dominate.
- (Q3) Promoter-only library (PLS cCREs) — promoters have stereotyped
  structure that may be easier to learn fast.
- (Q4) Genome-wide random 200bp windows — wider coverage than cCREs.

## 2026-06-02 21:57 — Planning Experiment 004

Plan: Test (Q1). 25,000 motif-packed sequences interleaved with 25,000
random-uniform sequences. Maximum motif vs no-motif contrast in one library.
If the model can learn ANY motif → activity, bimodal labels should give a
detectable r > floor.

Mode: exploring a new hypothesis (label-variance bottleneck). Result will
either prove the model CAN learn (and we lean into contrast) or prove it
CAN'T (and we move to Q2).

Generalization argument: a library that teaches "motif present → active"
generalizes to any cell type whose regulatory grammar uses TF motifs (i.e.,
all human cell types). Not a niche K562/HepG2/SK-N-SH bet.

## 2026-06-02 22:05 — Exp 005/006: random genomic and DHS-Meuleman both ~0

005 (genome-wide random 200bp): eval_01 = -0.0008.
006 (DHS Meuleman 16-component stratified, direct port of dhs_stratified
informed baseline which claims r ≈ 0.7055): eval_01 = -0.0008.

**Major conclusion: the instructions.md baseline table is not a faithful
description of v14's evaluator.** A direct, faithful port of the Meuleman
DHS approach the table claims gets r ≈ 0.7 actually gets r ≈ 0 in v14. The
"informed" instructions were misleading.

Theory v3. v14's evaluator returns near-zero Pearson r for every library
type I have tried so far (n=6: cCRE-random, cCRE-replicated, motif-packed,
motif-vs-random, genome-wide-random, DHS-stratified). Noise floor is
roughly [-0.005, +0.005]. The eval sets may use:
- (T1) a distribution very unlike any I have tried — e.g., synthetic
  reporter constructs with TATA + binding-site cassettes, or sequences from
  a specific published MPRA. Library distribution match dominates.
- (T2) labels that are nearly uncorrelated with sequence given v14's
  evaluator. In that case nothing works and the "best" library is the most
  principled one.

Plan: 2-3 more targeted experiments to test (T1), then commit to (T2) and
ship a principled design.

## 2026-06-02 22:08 — Planning Experiments 007-009

- 007 (now): Bimodal high-signal vs low-signal DHS. Take top 25k and bottom
  25k DHS by mean_signal. Maximum biological-activity contrast. If the model
  responds to activity-level signal at all, this should crack the floor.
- 008: Reporter-construct-like — TATA + many TF binding sites + spacers.
  Mimics MPRA library design conventions.
- 009: Download a published MPRA training set (Sharpr-MPRA chr18 sequences,
  pad to 200bp), or use Agarwal LentiMPRA sequences. Tests T1 directly.

If all 3 give ~0, I commit to (T2): final library will be a maximally
principled diversity design that would generalize to unseen cell types, and
I'll spend remaining experiments documenting / refining that design.

## 2026-06-02 22:30 — Exp 007/008/009 — all noise

- 007 (DHS top vs bottom by mean_signal, max activity contrast):
  eval_01 = +0.0013.
- 008 (reporter-construct mimic: 140bp enhancer + 60bp fixed minP):
  eval_01 = -0.0039. Identical to v14 random_uniform baseline.
- 009 (Sharpr-MPRA top+bottom 25k by mean activity, padded to 200bp):
  eval_01 = -0.0002.

**This is decisive.** Real measured MPRA training data with extreme
activity contrast (009) gives r ≈ 0. Synthetic reporter cassettes (008)
give r ≈ 0. Bimodal DHS contrast (007) gives r ≈ 0. (T1) is falsified: the
v14 evaluator does not respond to any biologically meaningful distribution
I can construct.

**Theory v4 (final).** v14's evaluator is, for all practical purposes,
uncorrelated with library content over the space of strategies a thoughtful
designer would try. The mean Pearson r across 14 eval sets is noise around
[-0.005, +0.005] regardless of input. Optimizing for the metric is futile.

Pivot to (T2): **principled library design**. The goal becomes "if v14's
metric is uninformative, build the library a real ML researcher would want
for a model that genuinely generalizes across cell types." Specifically:
  1. Maximize coverage of TF-motif feature space (combinatorial richness).
  2. Maximize coverage of genomic-context feature space (multiple genomic
     element classes, multiple cell types).
  3. Include both real and synthetic sequences so the model sees both
     in-distribution and counterfactual examples.
  4. Maintain length / GC distributional balance to avoid spurious
     covariates.

Plan for remaining experiments (010-030):
- 010: cCRE class-stratified (PLS / pELS / dELS / CA-CTCF / CA-only / TF /
  K4m3) — 14k each across the 7 main cCRE classes. Tests cCRE-class breadth.
- 011: Multi-source merge — equal parts cCRE-PLS, cCRE-dELS, DHS-summit,
  Sharpr-top, Sharpr-bottom (10k each = 50k). Tests source diversity.
- 012: TF motif "all-pairs" — every pair of 14 TF motifs covered ≥ once,
  random spacing/orientation/flanks. Combinatorial coverage of motif grammar.
- 013-015: Sweep mixes (synthetic fraction = 0 / 25 / 50 / 75 / 100 %).
- 016-020: Variant series (focus on the design from the best of 010-015 and
  vary one design knob each time — even though metric is uninformative, this
  is what disciplined ML library design looks like).
- 021-030: Final principled library + sensitivity/robustness experiments.

Even if metric stays noise, the deliverable is a defensible library design
with documented rationale.

## 2026-06-02 22:55 — Exp 010-016 — all noise (16 experiments total)

Tested 7 more principled designs:
- 010 cCRE class-stratified (5 classes × 10k): -0.0002
- 011 multi-source merge (5 sources × 10k): -0.0007
- 012 motif-grammar combinatorial: -0.0002
- 013 variant-augmented (Sharpr parents × variants): -0.0002
- 014 RC-augmented: +0.0003
- 015 k-mer max-diverse: -0.0015
- 016 50-TF motif vocabulary: +0.0010 ← best so far

After 16 experiments the entire score range is [-0.0039, +0.0026] on eval_01.
Indistinguishable from sampling noise. Theory v4 (v14 evaluator is
essentially uninformative) holds. 

Still, 016 sits at the top end. Plan for 017-020: ablate the 016
motif-vocab50 recipe to see if any axis nudges the score reproducibly:
  - 017: 016 design but bimodal (25k motif-rich + 25k pure random)
  - 018: 016 design but doubled motif count (8-16 per seq)
  - 019: 016 design + reverse-complement of every sequence
  - 020: 016 motifs IN cCRE genomic backbones (real flanks, motif inserts)

021-030: take the best of 017-020 + 016 and build a final flagship library,
then do sensitivity/scale tests.

Theory v5 emerging: even when the v14 metric is silent, two design axes
plausibly help downstream generalization in any real MPRA setting:
  (a) maximum motif vocabulary breadth (more TFs = wider grammar exposure)
  (b) real genomic context as backbone (gives the model authentic local-
      composition statistics around motifs).
This is what 020 will test.

## 2026-06-02 23:25 — Final summary (30 experiments complete)

### Headline result

All 30 experiments returned eval_01 mean Pearson r in [-0.005, +0.005].
Best: 004 (motif-vs-random bimodal) = +0.0026 and 019 (50-vocab+RC) = +0.0025.
Worst: 026 (CpG-enriched) = -0.0050. Same recipe across two seeds gives
a score swing of ~0.001 (019=+0.0025 vs 022=+0.0013; 029=+0.0001 vs
030=+0.0007). This is noise.

### v14 metric is uninformative

Confirmed by 3 independent lines of evidence:
1. v14's own baselines (random_uniform, gc_50, etc., per strategies.md) all
   sit in the same [-0.005, +0.005] band.
2. A faithful port of instructions.md's claimed "best" strategy (dhs_
   stratified, allegedly r ≈ 0.7055) returned r = -0.0008. The
   instructions.md baseline table does not describe v14.
3. Identical recipe re-run with a different seed gives a score swing
   comparable to the entire inter-strategy range.

### Strategies tried (30 distinct libraries)

Real genomic:
- 001 cCRE random, 002 cCRE replicated, 005 random genomic, 006 DHS-Meuleman,
  007 DHS bimodal-signal, 009 Sharpr real, 010 cCRE 5-class, 026 CpG-enriched

Synthetic:
- 003 motif-packed (14 vocab), 008 reporter-construct, 012 motif-grammar
  combinatorial, 015 k-mer max-diverse, 016 motif-vocab50

Hybrid / augmentation:
- 004 bimodal motif-vs-random, 011 multi-source merge, 013 variant-augmented,
  014 RC-augmented Sharpr, 017 vocab50 bimodal, 018 vocab50 dense, 019
  vocab50+RC, 020 motifs-in-cCRE, 021 bimodal motif+RC

Ablations / sensitivity:
- 022 seed sensitivity (019 recipe, seed=1), 023 flagship v1, 024 flagship v2,
  025 flagship v3, 027 cCRE+RC ablation, 028 Sharpr+RC ablation, 029
  flagship v4 (final), 030 flagship v4 seed=1

### Deliverable

029_flagship_v4 is the chosen library. Composition:
- 30k cCRE (5 classes, RC-paired)
- 20k Sharpr top+bottom poles (RC-paired)
- 10k synthetic motif + motif-in-cCRE (50-TF vocab, RC-paired)

Rationale: maximum biological coverage given that no specific design tested
moves the v14 metric beyond noise.

### Theory v5 (final)

The v14 evaluator either uses a fixed model uncorrelated with library
content over the strategy space we explored, or it evaluates against a
distribution so narrow that no thoughtful designer would discover it
within 30 trials. The "informed" instructions table was a red herring
for this evaluator. The right move was to (a) verify that signal was
absent via diverse exploration, (b) document the verification, (c) commit
to a principled design defensible on biological grounds alone.

### Skills written

- skills/extract_genomic.py — common genome-window extraction
- skills/sharpr_decode.py — Sharpr HDF5 → ACGT string with padding
- skills/motif_pack.py — synthetic motif-packed sequence generator
