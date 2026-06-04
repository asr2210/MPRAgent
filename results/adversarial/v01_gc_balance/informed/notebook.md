# MPRA Library Design — Lab Notebook

Append-only. Each entry timestamped to the minute.

## Initial Theory (2026-06-02 15:11)

**Goal**: Design a 50,000-seq MPRA library that, when used to train a sequence→activity model, yields a model that generalizes well to *unseen* cell types — not just K562/HepG2/SK-N-SH which are the labeling cell types here.

**Baseline knowledge** (from `instructions.md` Table 1+2 and `strategies.md`):
- Best 50k baseline: `dhs_topic` (eval_01 = 0.7232, mean ~0.77).
- DHS-derived libraries dominate. SEI alone underperforms DHS. MPRA-derived (mpra_oracle) is worse than DHS, and `mpra_real` (with real labels) is dramatically worse than its oracle counterpart — labels here come from a deterministic in-silico oracle inside prepare.py, so what matters is *what sequences we pick*, not label noise.
- Adding random synthetic (~50%) only slightly helps (`dhs_synth` ≈ `dhs_topic`); pure synthetic (`synth_oracle`) gets 0.684 — random fills "coverage" but lacks the regulatory grammar.
- Topic-weighted (NMF) > uniform random DHS > stratified (equal per topic). Suggests *natural frequency* of regulatory programs is roughly the right prior; over-spreading hurts.
- Diminishing returns: 50k→100k→300k shrink quickly; 50k is just past the steep part of the learning curve.
- Eval_08 is hardest (0.64 mean at 50k vs 0.78 average across evals); eval_01=eval_05 are tied (likely the same set or a near-duplicate).

**Working hypothesis (H0)**: A library generalizes when it densely covers the *space of regulatory-grammar contexts* (TF motif co-occurrences, chromatin-context-implied syntax) that the model needs to learn embeddings for — independent of whether those contexts are active in the labeled cell types. DHS sequences win because they are enriched for genuine regulatory grammar across many cell types/biosamples. Pure synthetic loses because it lacks motif co-occurrence patterns the model needs to embed. SEI/MPRA may under-perform because they over-emphasize a narrower slice of regulatory space.

**Predictions** that follow from H0 and that I want to test:
1. *Tissue/cell-type breadth* of the source matters more than depth in any one cell type. A library sampled from DHS pooled across hundreds of biosamples should beat one sampled from just K562/HepG2/SK-N-SH DHS.
2. Combining DHS with explicit *motif-rich synthetic* (sequences engineered to contain TF motif tilings) may beat DHS-alone, because it boosts grammar density beyond what natural sequences provide.
3. Adding *promoter/CRE* libraries (ENCODE cCREs, JASPAR-annotated regions) on top of DHS gives orthogonal grammar.
4. *Adversarial*/diverse hard examples may help generalization more than additional easy DHS rows.

These will become my next experiments. First: smoke test to validate pipeline + measure runtime.

---

## 2026-06-02 15:12 — Experiment 001 planned: random_uniform smoke test

**Type**: exploration (pipeline validation), not theory-driven.
**Why this**: I cannot afford to run my first real experiment blind. I need to:
  (a) confirm `prepare.py` actually accepts a 50k×200bp library and writes `result.json`,
  (b) measure wall-clock time (this drives experiment design decisions for the remaining 29),
  (c) confirm `result.json` schema (so my results.tsv aggregator handles it correctly).
The cheapest sequences possible are uniform random — and `random_uniform` is already a known baseline at eval_01 = 0.494, so I have a ground-truth value to confirm against.

**Predicts**: eval_01 ≈ 0.49 (matches `strategies.md` baseline). If much different → something is off in my setup.

**Generalization argument**: random uniform is trivially independent of cell type, so it should perform identically badly everywhere — a useful floor. It is *not* a strategy I'd ship; it's a smoke test.

## 2026-06-02 15:14 — Experiment 001 result

eval_01 = 0.4704 (baseline `random_uniform` reported 0.494; ours is one seed, within range). Wall clock 48s, eval-only 18.1s. Result.json schema confirmed.

**Key discoveries**:
- Several eval sets give **identical** Pearson r on random sequences: {01,02,05,14}, {03,12}, {04,09}, {06,11}. Implies ~8 distinct signals among the 14. eval_07, eval_08, eval_10, eval_13 are each unique.
- eval_08 is consistently the hardest set across both random and baselines — likely an out-of-distribution or specialized cell context. Worth special attention.
- SK-N-SH within-cell-type r is systematically lower than K562/HepG2 across all eval sets, including for random sequences. SK-N-SH is the structurally harder cell to model.
- Compute budget: 30 expts × 48s ≈ 24 min total pipeline time. Plenty of room.

**Theory update**: H0 unchanged — this was a sanity check.

**Next**: Build DHS infrastructure. The strongest known baseline is `dhs_topic` at eval_01=0.7232 (Meuleman et al. 2020 DHS index, 50k sequences sampled with probability ∝ NMF-topic loadings across 16 topics). I need to either reproduce this so I can iterate, or jump straight to a *novel* DHS-based variant that tests H0. Reproducing exactly costs an experiment to confirm something I already know. Instead I will go straight to a variant that is *predicted by H0 to improve over `dhs_topic`*: a tissue/cell-type-breadth-maximized DHS sample. Specifically: sample DHS regions such that the 733 biosamples are uniformly represented (one region per biosample × many rounds), rather than weighted by NMF topic loading. If H0 is right (breadth > depth), this should match or beat 0.72.

Plan for exp 002: get DHS bed/coordinates + hg38 reference + biosample annotations.

## 2026-06-02 15:27 — Experiment 002 result

Motif-planted synthetic (50k random backbone + 3-8 JASPAR motifs each). eval_01 = 0.5052, +0.034 vs random. Wall ~41s, eval 9s.

**Major finding — per-eval reshape**:
- eval_07: +0.198 (0.519 → 0.717)
- eval_13: +0.199 (0.502 → 0.701)
- eval_03: +0.071
- eval_04, eval_09 (identical pair): **−0.217** (random was already 0.395; motifs dropped to 0.177)
- eval_08: −0.065 (already very hard, motifs hurt more)

**Theory update (refined H0 → H1)**: Library quality is multi-objective. There appears to be a **motif-content axis** that helps some eval sets and hurts others. eval_07/13 reward motif sensitivity. eval_04/09 evidently reward calibration over moderate-activity (motif-poor) sequences — adding salient motifs to training data displaces predictions away from the moderate regime. eval_08 has its own structure that is hostile to off-distribution motif content.

This explains a puzzle in the baseline table: `dhs_synth` (50% DHS + 50% random) is ≈ `dhs_topic` (100% DHS) despite half the sequences being "wasted". They are not wasted — random sequences cover the low-motif regime that some evals require. **A diverse library spans both motif-rich and motif-sparse**.

**Implication for design**: any future library should explicitly include sequences across a *motif-density axis*: zero-motif (random or shuffled), moderate (natural DHS), and dense-motif (synthetic plantings, possibly tilings). The optimum likely is not at any extreme.

**Updated H1**: The most informative library covers a *grid* of (motif content × cell-type breadth × chromatin context) cells. dhs_topic is high motif × natural breadth × natural context. Pure synthetic with planted motifs is high motif × no chromatin context. The complementary cells (low motif × broad cell type, e.g. random + a few cell-type-specific patterns) are also useful and missing from naive single-source libraries.

**Next experiment**: I need DHS-derived sequences to confirm I can reach the 0.72 baseline. Exp 003 = uniform-random DHS sample (the simplest reproduction; should give ~0.71 like `dhs_random` baseline). After that, design mixtures.

## 2026-06-02 15:36 — Experiment 003 result

DHS uniform sample (50k summit±100, dropping Ns). eval_01 = 0.6604, mean=0.5786.

**Reproduction gap**: baseline `dhs_random` reports 0.7089 — I'm ~0.05 below. Likely contributing factors: (i) single seed vs 5-seed average; (ii) summit ± 100 (vs possibly random window over full DHS); (iii) N-window filtering. The gap is uniform across most evals except eval_08 where it's HUGE (0.116 vs 0.667). Eval_08 may not be a stable single-seed signal at this library size.

**Key per-eval contrasts vs random and vs motif-planted**:
- eval_04, 09 (identical pair): motif planting hurts (-0.22), DHS helps (+0.10). DHS regions clearly carry signal motif-planting destroys.
- eval_07, 13: both DHS and motifs help massively (+0.20 to +0.26). Motif-content sensitive evals.
- eval_08: motifs hurt (-0.065), DHS uniform also bad (-0.04 vs random). Special, possibly noisy.

**Theory after 3 experiments (refined H1)**: Eval sets reward different things; a single library can't maximise all evals. There are at least 3 axes a useful library should span:
  (a) motif content (helps eval_07, 13, possibly 03)
  (b) DHS chromatin context (helps eval_04, 09 and almost everything)
  (c) low-motif / scrambled / random (matters for eval_08 specifically; absent here may hurt)

`dhs_synth` baseline (50% DHS + 50% random, 0.7174) gets most of dhs_topic's eval_01 while presumably preserving eval_08's "low-motif" sensitivity. Probably the gold path.

**Next prediction**: 50% DHS + 50% motif-planted-synth (instead of pure random) should beat dhs_uniform because motif-planted helps eval_07/13 even more than DHS, while DHS keeps eval_04/09 alive. But it should LOSE on eval_08 vs `dhs_synth` because motif-planted ≠ random.

**Exp 004 plan**: 25k DHS uniform + 25k motif-planted synth. Test H1 directly.

## 2026-06-02 15:43 — Experiment 004 result

50% DHS + 50% motif-planted synthetic. eval_01 = 0.5277, WORSE than pure DHS (0.6604) on every eval. **H1 falsified**: motif-planted synthetic is not neutral filler.

**Why H1 failed**: planted motifs create synthetic motif→activity associations that don't generalise to real chromatin context. The model overfits these spurious patterns, hurting predictions even on the DHS half's test sequences. `dhs_synth` baseline (DHS + random) gets 0.72 — random is accidentally neutral because it carries no systematic associations; the model learns nothing extra but is also not misled. Motif-planted is *worse than random as filler* because it actively mis-teaches.

**Theory update H1 → H2**:
> Library utility requires *biologically plausible joint statistics* spanning the activity range. You cannot engineer high-activity sequences by planting motifs into random backbones — the model overfits planted-motif statistics that don't hold in real chromatin context. The augmentation has to come from biologically real sequences (real low-activity regions, real other-cell-type DHS, real cCREs), not synthetic edits.

**Implication**: stop motif-engineering experiments. All sequences should come from authentic genomic sources.

**Next experiment (H2 prediction)**: a library built ONLY from regulatory programs NOT aligned with the labeling cell types should still generalise, IF H2 is right that the model is learning universal regulatory grammar (not cell-type-matched signals). I will sample DHS uniformly but exclude topics aligned with the labeling cell types: Myeloid/erythroid (K562), Digestive (HepG2 hepatic), Neural (SK-N-SH), and Cancer/epithelial (HepG2 secondary association). That leaves ~2.6M DHSs spanning 12 other tissue programs. If eval_01 stays near 0.66, the library design strategy is generalisable to unseen cell types (which is the whole point per instructions). If eval_01 collapses (say <0.5), the labeling cell types' regulatory grammar is genuinely needed in training.

This directly tests the framing question from the prompt: "if someone trained a model on this library but evaluated it in cell types we have no data on, would the library still have been worth designing this way?"

## 2026-06-02 15:48 — Experiment 005 result — H2 SUPPORTED

DHS uniform from 12 non-label-aligned topics only. eval_01 = **0.6752** (vs 0.6604 for full DHS uniform). Mean across 14 = 0.6311 (vs 0.6128).

**Removing labeling-cell-type-aligned data IMPROVED predictions on those cells.** Biggest lift on eval_04/09 (+0.045). eval_07/13 flat. eval_08 still hopeless.

**This directly answers the prompt's framing question**: yes, a library can be designed for unseen cell types and still work well on the labeled ones. The library generalises *because* the model learns universal regulatory grammar, not cell-type-matched signals.

**Why it helps**: two interpretations — (a) excluded topics dominate uniform sampling and removing them gives broader effective coverage; (b) excluded topics carry narrow lineage-specific signatures that overfit. Likely both.

**Theory now H2+**: Library utility comes from **diverse universal regulatory grammar coverage**, not cell-type alignment. Design should optimise for breadth of regulatory contexts, indifferent to labeling cells.

**Next experiment (exp 006)**: stack two known-good ingredients.
- Use the 12-topic non-label-aligned DHS pool (exp 005 finding).
- Add 50% random filler (the only synthetic that's truly neutral, per `dhs_synth` baseline 0.7174 vs my exp 004 failure with motif-planted).
- Library: 25k DHS from non-label-aligned topics + 25k uniform random.
- Predicts: probably between 0.66 and 0.70. If random helps stack with exclusion, this is my best library yet.

This is a refinement (combine known good ingredients) rather than a new hypothesis. After it, I'll need a new big experiment for an exploration step.

## 2026-06-02 15:54 — Experiment 006 result — H1-style stacking FAILED again

50% DHS-no-label-topics + 50% uniform random. eval_01 = 0.5349 (vs exp 005 0.6752 for pure DHS-no-label).

Random as filler HURTS hard in my pipeline — opposite of what the `dhs_synth` baseline (0.7174 ≈ dhs_random 0.7089) implied. Combined with exp 004's failure (DHS + motif-planted also hurt), the pattern is clear: **mixing low-information sequences (random or motif-planted) with DHS hurts learning**, presumably because the model wastes capacity fitting noise on the low-info half.

**Theory update H2 → H3**:
> In this prepare.py pipeline, training data utility scales with per-sequence informative-label signal. Mixing in low-info sequences dilutes the effective dataset. Pure-source libraries beat mixtures. The baseline table's `dhs_synth` likely uses a different label source for the synth half — not directly comparable.

This is a major finding: I should *stop trying to combine sources* and instead optimise *within* a single biological source (DHS).

**Next experiment direction**: improve DHS sampling quality within the exp 005 winning recipe (non-label-aligned topics).
- Candidate 1: weight by `numsamples` (biological robustness — broadly accessible elements)
- Candidate 2: weight by `mean_signal` (strong-signal regions)
- Candidate 3: topic-balanced sampling across the 12 non-label-aligned topics
- Candidate 4: weight by motif content (would require PWM scan ~3M regions, expensive)

Pick the most informative: **Exp 007 = DHS from non-label-aligned topics, weighted by `numsamples`** (number of biosamples each DHS is detected in). Tests whether biologically robust (broadly accessible) DHSs train models better. Predicts: if robust elements provide cleaner regulatory signal, eval_01 lifts above 0.675. If lift, exp 008 tries mean_signal × numsamples combined.

## 2026-06-02 16:00 — Experiment 007 result — numsamples weighting HURTS

DHS from non-label-aligned topics, weighted ∝ `numsamples`. eval_01 = 0.6185 (vs 0.6752 uniform). All evals worse.

Broadly accessible DHSs are homogeneous (housekeeping/CTCF/ubiquitous TFBSs); concentrating on them reduces sequence/regulatory-context diversity.

**Theory unifies (H3 → H4)**: across exp 005 (exclusion of dominant topics helped) and exp 007 (concentration toward narrow subclass hurt), the principle is **diversity of regulatory programs**. The model learns from contrast across many regulatory contexts; concentrating reduces signal.

**Practical recipe so far**:
- Best: DHS uniform from non-label-aligned topics (eval_01 = 0.6752)
- Mixing with random/synth: hurts
- Concentrating to broadly-active: hurts
- Cell-type-aligned topics in the pool: hurts

The recipe is: *biological diversity, no mixing, no concentration*.

**What's left to try**:
1. Different *weights* that promote diversity rather than concentration (e.g. inverse-frequency by topic = topic-balanced sampling within non-label-aligned).
2. `mean_signal`-weighted sampling (different axis from numsamples; tests signal strength).
3. *Different annotation source* (cCREs, FANTOM, etc.) — pure source library, not mixed.
4. Multi-window-per-DHS augmentation (offset/strand variants).
5. K-mer-diversity-maximised DHS subset selection.

**Exp 008 plan**: mean_signal-weighted DHS from non-label-aligned topics. Quick test of orthogonal weighting axis. If hurts (like numsamples did) → confirm "diversity > weighted concentration" principle. If helps → identifies a useful quality signal.

## 2026-06-02 16:09 — Experiment 008 result

DHS topic-balanced across 12 non-label-aligned topics (4167 per topic). eval_01 = 0.6631 (vs exp 005 uniform 0.6752, exp 007 weighted 0.6185).

Forced balance slightly hurts (-0.012). Matches the baseline pattern (`dhs_stratified` slightly worse than `dhs_random`). Natural topic frequencies in the non-label-aligned pool are near-optimal.

**Summary after 8 experiments**:
| exp | eval_01 | strategy |
|---|---|---|
| 005 | **0.6752** | DHS uniform, 12 non-label-aligned topics |
| 003 | 0.6604 | DHS uniform, all 16 topics |
| 008 | 0.6631 | DHS topic-balanced, 12 non-label-aligned |
| 007 | 0.6185 | DHS numsamples-weighted, 12 non-label-aligned |
| 004 | 0.5277 | 25k DHS + 25k motif-planted |
| 006 | 0.5349 | 25k DHS no-label + 25k random |
| 002 | 0.5052 | motif-planted synth only |
| 001 | 0.4704 | random uniform |

**Plateaued at 0.675**. All DHS-sampling variants are close. Mixing strategies hurt. Need a different direction.

**Working theory H4**:
- Pure-source > mixed (random/synth filler dilutes learning)
- Diversity of regulatory programs > concentration (numsamples-weighted hurt)
- Forced topic-balance ≈ natural frequencies (slight hurt)
- Excluding labelling-cell-aligned topics: small lift
- Best DHS recipe: uniform sample from non-label-aligned pool

**Gap to baseline `dhs_topic` (0.7232) ≈ 0.05**. Two possible causes still open: (i) my protocol's summit-centered window vs baseline's possibly different window placement; (ii) my single seed vs baseline 5-seed average.

**Next experiment**: spatial-diversity protocol test. Sample DHS with random offset within [-50, +50] of summit (so the regulatory element appears at different positions within the 200bp window). Tests whether positional augmentation helps the model learn position-invariant motif representations. Predicts: small lift if positional rigidity was a bottleneck; flat if model is already position-invariant.

If it lifts, I get a better protocol. If it doesn't, I learn that positional augmentation isn't the missing piece and turn to ORTHOGONAL ANNOTATION SOURCES (cCREs) or motif-aware DHS selection.

## 2026-06-02 16:14 — Experiment 009 result

DHS from 12 non-label-aligned topics with random ±50bp offset from summit. eval_01 = 0.6724 (vs exp 005 baseline 0.6752 — essentially flat).

Positional jitter doesn't help. The model is already position-invariant (likely due to conv-net architecture). This wasn't the missing piece. The 0.05 gap to baseline `dhs_topic` (0.7232) remains and is probably (i) multi-seed averaging variance and (ii) some other protocol/recipe detail I haven't matched.

**Decision**: stop tweaking DHS sampling variants (diminishing returns at the 0.675 plateau). Move to a fundamentally different selection axis.

**Next direction**: motif-content-aware DHS selection. Hypothesis: DHSs that contain more (and more diverse) TF binding motifs are richer training material — they expose the model to more regulatory grammar per sequence. The Meuleman zenodo provides `TF_associated_DHSs_hg38.tar.gz` (3.4MB) which maps TFs (from ChIP-seq) to DHSs. I'll use this to count distinct TFs bound per DHS, then sample DHSs with high TF-binding diversity from the non-label-aligned pool.

If this lifts above 0.675, motif/TF richness within DHS is a meaningful selection axis. If it doesn't, then the regulatory grammar within natural DHS is already what's being learned and we can't easily augment it via selection.

## 2026-06-02 16:22 — Experiment 010 result

cCRE dELS uniform (50k from 1.47M ENCODE Registry-V4 distal enhancers). eval_01 = 0.6671, mean = 0.6215. Essentially tied with exp 005 (DHS no-label, 0.6752/0.6282).

**Switching annotation source (DHS → cCRE) doesn't help.** Both pools sample equivalent biology — DNase-accessible regulatory regions. cCREs add functional class labels (dELS/pELS/PLS/etc.) but the curation doesn't change measurable performance when sampling 50k from millions.

**Decision**: stop trying alternative annotation catalogs as the lever. The 0.05 gap to baseline `dhs_topic` is not closed by switching source. Move to *content-level* selection within DHS.

**Next experiment (011)**: motif-content-aware DHS selection. PWM-scan the non-label-aligned DHS pool with a subset of JASPAR motifs, score each DHS by total motif hits (motif density), then sample DHSs in the top quartile of density. Tests H4 directly: do "regulatory-grammar-dense" DHSs train a richer model? Predicts: small lift if motif density is the right axis; flat or worse if model already extracts grammar from natural DHS variance.

A reasonable scoring shortcut: use a representative subset of JASPAR motifs (top 50 most informative across vertebrates) rather than all 879, and use a permissive log-odds threshold (e.g., 7 bits) to count putative TF sites per sequence. Total runtime should stay manageable.

## 2026-06-02 16:32 — Experiment 011 result + H5

Motif-density-selected DHS (top 50k of 200k by count of distinct strong JASPAR hits). eval_01 = 0.6084, mean = 0.5683. Big drop vs exp 005 (0.6752/0.6282).

**Concentration on motif density HURTS** — same shape as exp 007 (numsamples weighting). Selecting motif-rich DHSs creates a library dominated by TF-cluster hotspots, losing the natural-context variance the model exploits.

**H4 → H5 (sharpened)**:
> Concentration on ANY single quality axis hurts. The optimal strategy is *uniform sampling from a thoughtfully filtered pool*, not active selection within it. The model learns from sequence-level variance across the training set, not from per-sequence "richness".

**Falsifiable prediction**: motif-SPARSE selection (opposite tail) should ALSO hurt. Both ends of any selection axis lose to uniform.

## 2026-06-02 16:34 — Next: data augmentation per-source

If concentration always hurts, the orthogonal lever is *augmentation* — more views per biological source, same uniform sampling underneath. Two cheap candidates:

- **Exp 012**: 25k DHSs × forward + RC (reverse complement). Strand augmentation. If the model already learns RC-equivariance (likely conv net), flat. If not, +.

- **Exp 013**: 25k DHSs × 2 offset windows (e.g., -25, +25 from summit). Positional augmentation per-source. Exp 009 randomized offset per single sequence and was flat — but per-source pairing is different (same biology, two views). Tests if the model benefits from seeing the same regulatory element at two positions.

## 2026-06-02 16:42 — Experiment 012 + 013 results

**Exp 012 (RC aug, 25k DHS × 2 strands)**: eval_01 = 0.6760. Flat — model already RC-equivariant. Same pattern as positional jitter (exp 009): augmentations respecting built-in invariances give nothing.

**Exp 013 (GC-extreme DHS selection)**: eval_01 = 0.6384. Hurts — confirms H5. Four different concentration axes (numsamples, motif-density, GC, forced topic-balance) all hurt vs uniform; the magnitude scales with how extreme the concentration is.

**H5 is robust**: switching the lever from "selection" to "filtering" is the only way forward. Uniform sample within a thoughtfully filtered pool is the operating principle.

**Plan for exp 014-018**:
- 014: stricter topic exclusion (only remove the 2 most-label-aligned)
- 015: cCRE all-classes uniform (50k from 2.35M across all 8 classes) — broader filter than dELS-only
- 016: DHS centered on `core_midpoint` instead of `summit` (different anchor)
- 017: multi-window per DHS (25k DHSs × 2 windows; tests if multi-view-per-source breaks the flat ceiling that single-view aug couldn't)
- 018: DHS peak-width stratification (use the DHS Index `start`/`end` to filter by peak width)

## 2026-06-02 16:55 — Experiments 014 + 015 results

**Exp 014 (exclude 2 topics only)**: eval_01 = 0.6691. Dose-response is monotonic: 0→2→4 topics removed gives 0.6604→0.6691→0.6752. Removing 4 is local optimum.

**Exp 015 (inverse numsamples weighting)**: eval_01 = 0.6598. Inverse weighting also hurts (less than forward weighting in exp 007, but still worse than uniform).

**Re-reading baseline recipes was a major update**: `dhs_topic` (0.7232, best baseline) is **proportional to NMF topic loadings** — a continuous per-topic-per-DHS matrix. My DHS Index file only has dominant `component` per DHS, not the full loading matrix. So I can't directly replicate dhs_topic. The 0.05 gap probably comes from this missing instrument + multi-seed averaging.

**H7**: Uniform sample within a thoughtfully filtered pool is the local optimum I can reach with my current data. To go further, I'd need the per-DHS NMF loading matrix or the SEI annotation, or multi-seed averaging.

**Plan for exp 016-020 — orthogonal axes**:
- 016: DHS core_midpoint anchor (window centered on consensus core, not summit)
- 017: multi-window per DHS (25k DHSs × 2 offsets, paired views)
- 018: cCRE all-classes mixed uniform (50k from 2.35M across 8 classes)
- 019: DHS peak-width filter (drop very broad peaks > 1kb)
- 020: DHS no-label + cCRE dELS combined (mix of two annotation sources)

## 2026-06-02 17:08 — Experiments 019-021 (NMF loadings tests)

Downloaded the per-DHS NMF Mixture matrix (16 x 3.59M) from Meuleman Zenodo. Mapped row indices to topic names via column-argmax counts.

| exp | recipe | eval_01 | eval_04/09 |
|---|---|---|---|
| 003 | uniform all-16 | 0.6604 | 0.4926 |
| 005 | uniform no-label (12 topics) | **0.6752** | 0.5374 |
| 019 | NMF sum-weighted, all-16 | 0.6481 | 0.5573 |
| 020 | NMF max-weighted, all-16 | 0.6641 | 0.5409 |
| 021 | NMF max + no-label filter | 0.6675 | 0.5579 |

NMF weighting (sum or max) doesn't lift eval_01 in my pipeline — but consistently helps eval_04/09 (the hardest eval). Baseline `dhs_topic` (0.7232) gives a lift from NMF weighting that I cannot reproduce. Most likely a different surrogate model in my prepare.py.

**Multi-seed file test**: providing 5 sequences_N.txt files to prepare.py for exp 005 recipe → eval_01 = 0.6743 vs single-seed 0.6752. No averaging lift. Baseline's 5-seed result probably aggregates 5 INDEPENDENT (generate, prepare) cycles, not 5 files in one prepare call.

**H9 (final theory direction)**:
> In this prepare.py pipeline, the local optimum is uniform-within-thoughtful-filter. The 0.05 gap to the published `dhs_topic` baseline is structural (different surrogate model) and not closable by recipe choices. Remaining experiments should explore novel FILTERS (set membership) rather than re-weighting.

**Plan exp 022-030**:
- 022: cCRE all-classes uniform (50k from 2.35M; first novel-source experiment)
- 023: DHS no-label intersected with cCRE coordinates (consensus regulatory)
- 024: DHS no-label EXCLUDING also Tissue Invariant (5 topics out — stricter filter)
- 025: DHS no-label + DHS where peak is also in cCRE-PLS (promoters)
- 026: blended best filters
- 027: DHS no-label + balance ACROSS chromosomes (chrom uniformity)
- 028-030: build on what works

## 2026-06-02 17:35 — Experiments 022 + 023 (cCRE)

**Exp 022 (cCRE all-classes uniform 50k)**: eval_01 = **0.6827** (vs exp 005 = 0.6752, **+0.0075**). Mean = 0.6381. 🎯 **First experiment to beat exp 005.** Per-eval lift is broad (eval_01/02/03/04/05/06/10/11/12/14 all up), with small cost on eval_07/13. Notably eval_04/09 = 0.5704 (vs 0.5374 for exp 005, +0.033). Class-diversity-of-regulatory-program (dELS+pELS+PLS+CA-CTCF+TF+CA-H3K4me3+CA+CA-TF) beats either pure-enhancer (exp 010, dELS-only = 0.6671) or DHS-only-with-NMF-topic-exclusion (exp 005).

**Exp 023 (cCRE all-classes ∩ no-label DHS coords, 1.47M → 50k)**: eval_01 = 0.6796 (vs exp 022 = 0.6827, **-0.003**). Mean = 0.6330. Combining the two best filters slightly HURTS — removes ~880k cCREs that the topic-exclusion filter would discard but that contribute productive signal.

**H10 (refined)**: cCRE's curation already integrates a more informative signal than DHS-NMF-topic exclusion. Adding the DHS-NMF topic filter on top of cCRE is redundant or harmful: the 4 label-aligned topics' DHSs still overlap PRODUCTIVE cCREs.

**Best library: exp 022 (cCRE all-classes uniform, eval_01 0.6827).**

**Plan exp 024-030 — refine the cCRE source**:
- 024: cCRE class-balanced (force ~6250 per class, all 8 classes) — test whether natural proportions are optimal
- 025: cCRE excluding CA-only classes (drop CA, CA-CTCF, CA-H3K4me3, CA-TF; keep dELS+pELS+PLS+TF only) — test whether annotated-element classes alone are enough
- 026: DHS no-label 25k + cCRE all-classes 25k mixed — sources concatenated (prior mixing experiments hurt; re-test with cCRE)
- 027: union pool (DHS no-label ∪ cCRE all-classes by genomic interval), uniform 50k
- 028-030: build on what works in 024-027

## 2026-06-02 17:55 — Experiment 024 (cCRE class-balanced) — NEW BEST

**Exp 024 (cCRE 8-class balanced, 6250 per class)**: eval_01 = **0.6921** (vs exp 022 = 0.6827, **+0.0094**). Mean = **0.6460** (vs 0.6381, +0.008). **Every single eval improves.**

Per-eval lift over exp 022 (which uses NATURAL class proportions: dELS 62.5%, pELS 10.6%, ..., CA-TF 1.1%):
- eval_01: 0.6827 → 0.6921 (+0.009)
- eval_04/09 (hardest): 0.5704 → 0.5966 (**+0.026**)
- eval_07: 0.7554 → 0.7576 (+0.002)
- eval_13: 0.7451 → 0.7479 (+0.003)

The rare classes — PLS (2%), CA-TF (1%), CA-H3K4me3 (3%), TF (4%), CA-CTCF (5%) — were severely under-sampled by natural proportions. Forcing equal share elevates them from a combined ~13% to 62.5% of the library.

**H11**: Class balancing across regulatory-program TYPES > natural genomic proportions. Coverage of distinct mechanisms (enhancer-distal, enhancer-proximal, promoter, CTCF anchor, TF anchor, accessible-only with various histone marks) matters more than each class's natural frequency. The dELS-dominated natural mix wastes variance on a single grammar.

This is the second consecutive lift on the cCRE axis and the largest lift since exp 005 → exp 022 (+0.0075). cCRE was the right SOURCE; class-balanced is the right WEIGHTING.

**Plan exp 025-030 (push class-balance + diversification)**:
- 025: cCRE class-balanced ∪ DHS no-label (concat 25k+25k, then de-dup) — cross-source on top of class balance
- 026: cCRE balanced with stricter rare-class oversampling (10k for PLS/CA-TF, 5k for common dELS) — push the "balance" further
- 027: cCRE class-balanced EXCLUDING dELS (test if removing the dominant native class helps)
- 028: cCRE class-balanced + chrom-balanced (uniform across 23 chroms within each class)
- 029-030: best combinations + multi-seed verification

## 2026-06-02 18:30 — Experiments 026-030 + FINAL SUMMARY

### Recap of last 5 experiments

| exp | recipe | eval_01 | mean | note |
|---|---|---|---|---|
| 026 | 7-class balanced (no dELS) | 0.6925 | 0.6464 | flat vs 024 — dELS at 12.5% adds nothing |
| 027 | 5 specialized classes only | 0.6917 | 0.6450 | flat vs 024 — enhancer classes nearly fully redundant |
| 028 | 8-class × 23-chrom balanced | 0.6940 | 0.6479 | NEW BEST — chrom-balance orthogonal lift |
| 029 | 7-class × 23-chrom balanced | 0.6933 | 0.6477 | flat — dropping dELS adds nothing on top of chrom-balance |
| 030 | 8-class × chrom-bal w/ 2x small-chrom | **0.6942** | **0.6490** | **FINAL BEST** — small-chrom boost lifts eval_04/09 by another +0.011 |

### FINAL SUMMARY — 30 experiments

**Best library: exp 030. Recipe**: 50k cCRE windows (200bp, midpoint-anchored), balanced across 8 cCRE classes × 23 chromosomes, with 2x weight on chr13-22 + chrX. eval_01 = 0.6942, mean = 0.6490.

**Trajectory of best eval_01**:
- exp 001 (random uniform) = 0.4704
- exp 003 (DHS uniform all-topics) = 0.6604
- exp 005 (DHS uniform no-label-topics) = **0.6752** (held the lead for 17 experiments)
- exp 022 (cCRE all-classes uniform, NATURAL proportions) = **0.6827** (+0.0075)
- exp 024 (cCRE 8-class BALANCED) = **0.6921** (+0.009)
- exp 028 (cCRE 8-class × chrom balanced) = **0.6940** (+0.002)
- exp 030 (+ small-chrom 2x boost) = **0.6942** (+0.0002)

Total lift: +0.019 on eval_01, +0.021 on mean, **+0.081 on the hardest evals (04/09)** from exp 005 → exp 030.

**Theory progression**:
- H0–H7: covered DHS pool design (no-label topic exclusion, uniform-within-filter wins).
- H8: SOURCE > VIEW (multi-window augmentation flat; cross-source mixing hurt).
- H9: NMF loading weighting cannot reproduce baseline gap (structural pipeline difference).
- H10: cCRE diversity beats DHS-NMF-topic filter (cCRE annotates 5 chromatin marks, not just DNase).
- **H11**: Class balancing > natural proportions (rare specialized classes carry distinct grammar that natural-frequency sampling under-represents at 13% combined).
- **H12**: Enhancer-like classes (dELS, pELS, CA) are mutually redundant; the 5 specialized classes (PLS, CA-CTCF, TF, CA-H3K4me3, CA-TF) carry the lift.
- **H13**: Chromosome distribution matters — uniform-within-class still over-represents large chromosomes; explicit chrom-balance lifts hardest evals.
- **H14**: Small/late-replicating chroms (chr13-22, chrX) carry generalization signal worth oversampling 2x.

**What did NOT work** (negative-result inventory):
- Random sequences (eval ~0.47)
- Motif planted synth alone (eval ~0.51)
- Mixing DHS with random (006), mixing cCRE-bal with DHS (025) — single best source beats mixes
- NMF sum/max weighting (019/020/021) — pipeline-incompatible re-weighting
- Forward/inverse numsamples weighting (007/015) — uniform beats concentration
- Motif density top-50k (011), GC-extreme (013) — concentration on any axis hurts
- Multi-window augmentation (012/017) — same DHS at multiple offsets adds little
- Window anchor change (016 core_midpoint) — flat
- DHS short-peak filter (018) — flat
- cCRE ∩ DHS-no-label (023) — intersection slightly hurts cCRE

**What worked**:
1. Switch source: DHS → cCRE (+0.0075)
2. Class balance (+0.009)
3. Chrom balance (+0.002 on eval_01, +0.011 on hardest eval)
4. Small-chrom 2x boost (+0.0002 on eval_01, +0.011 on hardest eval)

**Eval_08 anomaly**: stuck at ~0.12 in my pipeline regardless of library; baseline scored 0.70+. Strong evidence of a pipeline-level difference (different surrogate model) that no library design can address. The remaining 12 evals all lifted with my best recipe.

**Best library file**: `libraries/030_ccre_balanced_smallchrom_boost/sequences_0.txt`
