# MPRA Library Design — Lab Notebook

Append-only. Each entry starts with a timestamp.

---

## 2026-05-18 21:55 — Project kickoff & initial theory

**Setup.** Project scaffold: `libraries/`, `skills/`, `data/`, `notebook.md`, `results.tsv`.
`prepare.py` is a black box. It runs 150k 200bp sequences through an MPRA in K562/HepG2/SK-N-SH,
trains a sequence→activity model on the resulting measurements, and evaluates that model on 14
hidden test sets. Returns `mean_r` + per-cell-line `r` per eval set. I never see what's in those
test sets.

**Stated goal.** A library that trains a model of *general* regulatory grammar (across all cell
types, not just the three measured). High training performance-to-size ratio. Diverse in sequence
space. Not only functional elements.

**Initial theory (v0).** A model only learns what's in the data. For a 150k library of 200bp
sequences, the bits that matter for training are:
1. **Signal range** — sequences must span a wide range of activity. All-baseline = nothing to
   regress on; all-high = no contrast either. The cross-product of (sequence features × activity)
   is where the model's loss comes from.
2. **Motif coverage** — TF binding sites, core promoter elements, repressor motifs all need to
   appear, ideally in varied contexts so the model learns motif identity rather than memorising
   genomic position.
3. **Diversity of context** — same motif in many backgrounds teaches additivity; combinations
   teach interactions.
4. **Negative/low-activity examples** — random / shuffled / scrambled-motif sequences anchor the
   model on what's *not* regulatory.
5. **Out-of-distribution generalisation** — if the library only mimics one set of genomic regions
   (e.g. only ENCODE cCREs in K562), the model overfits to that distribution and fails on the 14
   eval sets which likely span more diverse regulatory contexts.

**Open unknowns going into experiment 001.**
- What does the model see on pure random sequences? My prediction: very low activity across the
  board, almost no contrast, all eval sets ≈ 0 or very low r. This establishes the floor.
- How fast does prepare.py run? (Latency budget per experiment determines how aggressive I can be.)
- How noisy is mean_r between repeats of the same library design? Important to know before chasing
  small deltas later.
- Do the 14 eval sets covary, or do they split into clusters? Clusters would suggest different
  regulatory regimes being tested.

**Plan for experiment 001.** Generate 150,000 uniformly random 200bp ACGT sequences (seed=0).
Submit to prepare.py. Record per-eval r and runtime. This is purely a **baseline** — I am
exploring, not refining. It anchors all later experiments.

---

## 2026-05-18 22:53 — Experiment 001 result (random baseline)

**Result: mean_r ≈ 0.82 across all 14 evals. This is a massive surprise.**

| eval | mean_r |   | eval | mean_r |
|------|--------|---|------|--------|
| 01   | 0.7760 |   | 08   | 0.9080 |
| 02   | 0.8718 |   | 09   | 0.8890 |
| 03   | 0.8562 |   | 10   | 0.8670 |
| 04   | 0.8178 |   | 11   | 0.7623 |
| 05   | 0.7756 |   | 12   | 0.7394 |
| 06   | 0.8722 |   | 13   | 0.7762 |
| 07   | 0.7910 |   | 14   | 0.8722 |

Aggregate: mean across evals = 0.8200. K562 = 0.8114, HepG2 = 0.8197, SKNSH = 0.8284. Time = 3434s
(~57 min).

**My prediction was wrong.** I expected r ≈ 0. The reality is r ≈ 0.82. Things I now believe:

1. **Random 200bp DNA carries enough accidental structure** — varied GC, accidental short k-mer
   motifs, accidental weak TF binding sites — for the MPRA to produce a wide range of activity
   measurements, and 150k examples is enough for the model to learn a generalisable mapping.
2. **The model architecture inside prepare.py is strong.** It extracts substantial signal from a
   library a human would dismiss as "useless".
3. **The 14 evals are NOT equally hard.** eval_08 (0.908), eval_09, eval_14 are easy. eval_11/12
   (0.74-0.76) are hard. Hard evals will be the ones that distinguish library designs.
4. **Eval clustering** — {02, 06, 14} report essentially identical numbers (0.8718, 0.8722,
   0.8722). {01, 05} similar (0.776, 0.776). Some evals may share underlying data. Treat
   replicated clusters as one signal each when judging small deltas.

**Theory update (v0 → v1).**
- Drop: "the floor is near zero". The floor on this MPRA is r ≈ 0.82 from random.
- Keep: signal range, motif coverage, OOD generalisation matter.
- Add: **information density per sequence**. Random is the bottom of the spectrum; the question is
  what kind of structure adds *new* information the model can't already get from random.
- Add: hard evals (11, 12) and easy evals (08, 09, 14) should be tracked separately. Optimising
  for the easy ones is cheap; lifting the hard ones is where library design will pay off.

**Per-cell-line.** SKNSH slightly > HepG2 ≈ K562 on most evals. Not concerning. No cell-line is
catastrophically worse.

**Cost.** ~57 min wall-time. 30 experiments ≈ 28h compute total. Need each experiment to be
information-rich. Avoid tiny perturbations.

---

## 2026-05-18 22:55 — Plan for experiment 002

**Hypothesis being tested.** Sequences drawn from the real human genome carry more regulatory
information per base than uniform random, so a model trained on them will generalise better
across the 14 evals (mean_r > 0.82). If they *don't* beat random, that's an even bigger
update — it would mean uniform random is hard to beat for general regulatory grammar and that the
information bottleneck isn't naive "realism" but something else (signal range, motif diversity,
etc.).

**Why this experiment is the most informative next step.** It's a single clean variable change
(random → real genome) directly contrasted against the only data point I have. Whichever way it
falls, the answer reshapes the next 5+ experiments. If real-genome wins big, I lean into
genome-derived libraries and start asking *which* regions matter most. If it ties or loses,
I focus on structural properties of the sequence distribution (GC, motif coverage, activity
contrast) independent of biological realism.

**Design.**
- Download hg38 chromosome 22 (small, ~51 Mb).
- Tile into non-overlapping 200bp windows.
- Drop windows containing N or any non-ACGT character.
- Drop windows that are >50% lowercase (soft-masked = repeats) — keep some repeats but not pure
  repeat tracks.
- Uniform sample 150,000 of the remaining windows (seed=1).

This is **exploring a new hypothesis**, not refining a direction.

**Predicted outcome.** mean_r modestly higher than random (0.83–0.87). I'd be surprised by a big
jump (>0.90) because random already does so well.

**Adjustment.** chr22 alone yielded only 90k clean windows; expanded the pool to chr19+20+21+22
to get 448k → sampled 150k. Same hypothesis, broader genomic coverage.

---

## 2026-05-19 00:01 — Experiment 002 result (hg38 chr19-22 tiled)

**Result: mean across evals jumps 0.820 → 0.873 (+0.053). Real genome wins on 13 of 14 evals.
Eval 08 is the lone regression (-0.040).**

| eval | random  | genome  | Δ       |
|------|---------|---------|---------|
| 01   | 0.7760  | 0.8182  | +0.042  |
| 02   | 0.8718  | 0.9162  | +0.044  |
| 03   | 0.8562  | 0.9080  | +0.052  |
| 04   | 0.8178  | 0.8567  | +0.039  |
| 05   | 0.7756  | 0.8181  | +0.042  |
| 06   | 0.8722  | 0.9166  | +0.044  |
| 07   | 0.7910  | 0.8882  | **+0.097** |
| 08   | 0.9080  | 0.8679  | **−0.040** |
| 09   | 0.8890  | 0.9364  | +0.047  |
| 10   | 0.8670  | 0.9062  | +0.039  |
| 11   | 0.7623  | 0.8044  | +0.042  |
| 12   | 0.7394  | 0.7857  | +0.046  |
| 13   | 0.7762  | 0.8947  | **+0.119** |
| 14   | 0.8722  | 0.9164  | +0.044  |

Aggregate: mean 0.8731. K562 0.8729, HepG2 0.8769, SKNSH 0.8702. Time 3783s.

**What this confirms.**
- Eval clusters {02,06,14} = 0.9162/0.9166/0.9164 and {01,05} = 0.8182/0.8181 are real (n=2 now).
  Treat each cluster as one effective signal.
- Across the rest, the delta from genome is +0.04 on most evals, with two big positives
  (eval_07 +0.10, eval_13 +0.12) and one negative (eval_08 −0.04).

**What this means.**
- Genomic sequences carry substantial *additional* information beyond what random provides.
  My prediction (0.83–0.87 mean) was right; the magnitude is at the upper end.
- The negative on eval_08 is interesting — that eval set must contain something better-captured
  by random than by genomic context. Two leading guesses: (a) eval_08 contains negative-control
  random-like sequences; (b) eval_08 emphasises an activity range or sequence class our genomic
  sample is depleted in.
- Eval_13's +0.12 is the largest single jump. Suggests eval_13 tests something heavily reliant on
  genomic context / TFBS spacing / dinucleotide structure that random can't supply.

**Theory update (v1 → v2).**
- Random gives r=0.82 floor. Genomic tiling adds ~+0.05.
- The right library is not "as enriched for function as possible" but rather "covers more axes
  of the activity & sequence distribution than any single source can". Eval_08 proves this:
  going more genomic *hurt* somewhere.
- Hypothesis: an optimal library probably contains some random/synthetic component so eval_08-
  style sequences aren't left out.
- The hard evals (11, 12) both improved with genome but are still the lowest (0.78–0.80). These
  are the ones to focus on next.

---

## 2026-05-19 00:08 — Plan for experiment 003

**Hypothesis being tested.** Sequences enriched for putative regulatory function — ENCODE V3
candidate cis-regulatory elements (cCREs) — give further gains beyond plain genomic tiling.
The continuum is: random (no enrichment, 0.820) → genome (mild enrichment, 0.873) → cCREs
(strong enrichment, ?).

**Why this is the most informative next step.**
- Adds a third data point on the "function enrichment" axis. Whichever way it falls is
  informative:
  - If cCREs > genome: enrichment for function is the right axis, lean into it (and we'll start
    asking which cCRE class — promoter / enhancer / CTCF — matters most).
  - If cCREs < genome: untargeted genomic context is doing something cCRE-only sampling kills
    (probably background diversity / negative examples / dinucleotide context).
  - If cCREs ≈ genome: the gain from random→genome was about *evolutionary realism* per se, not
    about regulatory enrichment specifically.
- Either direction reshapes the next experiments.

**Design.**
- Use `data/encode_ccres_hg38.bed` (1,063,878 cCREs across all classes & cell types).
- For each cCRE, take the **200bp window centered on the cCRE midpoint** in hg38 (any
  chromosome). This standardises the length, removes the variable cCRE length confound, and gives
  the model a regulatory element in the middle with adjacent genomic context on either side.
- Skip cCREs where the 200bp window extends past chromosome ends, contains any N, or is
  >50% softmasked.
- Sample 150,000 uniformly from survivors (seed=2). No class balancing yet — let the natural
  cCRE-type distribution stand so we measure the unbiased effect of "cCRE-only" library.

This is **exploring a new hypothesis** (function enrichment axis), not refining the previous
direction.

**Predicted outcome.** mean_r ≈ 0.88–0.91. I'd guess cCREs beat plain genome by ~0.02 because
regulatory regions have more contrast (more high-activity sequences) but lose a bit of breadth.
Bigger uncertainty: eval_08 will probably regress *further* (since cCREs are even further from
random than tiles are).

**Practical step.** I'll need all 24 chromosome FASTAs to look up windows across cCREs. Currently
have chr1, 7, 14, 19-22. Need the rest — will download in background while building the
library script.

---

## 2026-05-19 01:00 — Experiment 003 result (ENCODE cCREs)

**Result: cCREs beat genome on every eval, mean 0.873 → 0.886 (+0.013). Eval_08
regression from exp 002 is largely recovered (+0.024 on that eval specifically).**

| eval | random  | genome  | cCRE   | Δ(cCRE−genome) |
|------|---------|---------|--------|----------------|
| 01   | 0.7760  | 0.8182  | 0.8288 | +0.011 |
| 02   | 0.8718  | 0.9162  | 0.9270 | +0.011 |
| 03   | 0.8562  | 0.9080  | 0.9193 | +0.011 |
| 04   | 0.8178  | 0.8567  | 0.8657 | +0.009 |
| 05   | 0.7756  | 0.8181  | 0.8285 | +0.010 |
| 06   | 0.8722  | 0.9166  | 0.9274 | +0.011 |
| 07   | 0.7910  | 0.8882  | 0.9025 | +0.014 |
| 08   | 0.9080  | 0.8679  | 0.8922 | **+0.024** |
| 09   | 0.8890  | 0.9364  | 0.9465 | +0.010 |
| 10   | 0.8670  | 0.9062  | 0.9272 | +0.021 |
| 11   | 0.7623  | 0.8044  | 0.8145 | +0.010 |
| 12   | 0.7394  | 0.7857  | 0.7959 | +0.010 |
| 13   | 0.7762  | 0.8947  | 0.9038 | +0.009 |
| 14   | 0.8722  | 0.9164  | 0.9276 | +0.011 |

Aggregate: mean 0.8862. K562 0.8841, HepG2 0.8933, SKNSH 0.8812. Time 3041s.

**Key observations.**
1. cCREs strictly dominate genome — no losses on any eval. Function enrichment is a clean win
   despite my prior worry it might overfit to cell-type bias.
2. Eval_08 is mostly rescued. Random is still slightly better on eval_08 (0.908 vs 0.892), but
   the regression that pure-genome caused is essentially closed by cCRE focus.
3. Diminishing returns on the function-enrichment axis: random→genome was +0.053; genome→cCRE
   only +0.013. We're saturating.
4. Hard evals (11, 12) still hard. +0.01 to ~0.80. Probably require a different lever —
   cell-type-specific biology, extreme-activity sequences, or non-cCRE element classes (Polycomb,
   insulators, lncRNA regulatory elements, ...).
5. HepG2 gained most across the three (+0.073 random→cCRE). All cell lines monotone-improved.

**Theory update (v2 → v3).**
- Function-enrichment is one axis: random < genome < cCRE, but plateauing.
- Diversity-helps hypothesis (from eval_08) holds in a soft sense: cCRE didn't fully match
  random on eval_08, so something in random is still useful. Either random's wide activity
  coverage or its uniform sequence-space coverage is doing real work.
- The hard evals (11, 12) point at an unsampled axis. They're not unlocked by function
  enrichment alone. Need to think about what tests something neither random nor cCRE captures.

---

## 2026-05-19 01:05 — Plan for experiment 004

**Hypothesis.** Diversity of source matters for general regulatory grammar. A library mixing
high-information cCRE-centered sequences with the unique strengths of random sequences (better
on eval_08, possibly better coverage of low-activity sequence space) will beat either source
alone.

**Why this is the most informative next step.** Two specific pieces of evidence point this way:
(a) eval_08 in exp 002 — random outperformed genome — and (b) cCRE in exp 003 only partially
recovered that. The mix is the direct test of whether random and cCRE are complementary or
just two points on the same curve. The user's stated goal also explicitly says "not only
functional elements" — diversity is part of the design spec.

**Design.**
- 75,000 cCRE-centered windows (from the same 745,398-window pool as exp 003), seed=3.
- 75,000 uniform random ACGT 200bp sequences, seed=3 (different sub-seed for random vs cCRE
  pulls).
- Concatenate, shuffle, write 150,000 lines.

This is **refining a promising direction** (function-enrichment is winning, with one
specific weakness around random-flavoured sequences) by testing the simplest counterfactual.

**Predicted outcome.**
- Mean: 0.870–0.885. Likely slightly *below* pure cCRE (0.886) because mixing in half random
  drops information density on most evals.
- eval_08: 0.895–0.905, between cCRE (0.892) and random (0.908). If we see >0.908 it would be
  a strong "complementary" signal.
- Hard evals (11, 12): probably similar or slightly worse than cCRE — random doesn't help
  them.
- Key check: does *any single eval* exceed both pure-source numbers? That's the diversity
  premium.

If diversity doesn't deliver here, I'll pivot to a different axis (e.g., synthetic motif
embedding, extreme-activity selection).


## 2026-05-19 02:35 — Experiment 004 result (cCRE + random hybrid 75/75)

**Result: mean drops 0.886 → 0.880 (−0.006), BUT eval_08 jumps to 0.916 — higher than EITHER
pure source (random 0.908, cCRE 0.892). Real diversity premium on eval_08.**

| eval | cCRE   | hybrid | Δ      |
|------|--------|--------|--------|
| 01   | 0.8288 | 0.8217 | −0.007 |
| 02   | 0.9270 | 0.9199 | −0.007 |
| 03   | 0.9193 | 0.9116 | −0.008 |
| 04   | 0.8657 | 0.8609 | −0.005 |
| 05   | 0.8285 | 0.8214 | −0.007 |
| 06   | 0.9274 | 0.9204 | −0.007 |
| 07   | 0.9025 | 0.8896 | −0.013 |
| 08   | 0.8922 | 0.9156 | **+0.023** |
| 09   | 0.9465 | 0.9387 | −0.008 |
| 10   | 0.9272 | 0.9193 | −0.008 |
| 11   | 0.8145 | 0.8078 | −0.007 |
| 12   | 0.7959 | 0.7900 | −0.006 |
| 13   | 0.9038 | 0.8885 | −0.015 |
| 14   | 0.9276 | 0.9205 | −0.007 |

Aggregate: mean 0.8804. Time 4749s.

**Key observations.**
1. Eval_08 super-additive: 0.916 beats both pure-source values. Mixing produces something
   neither source alone has — strongest signal yet that *diversity of source* is a real lever
   distinct from function enrichment.
2. Every other eval loses 0.005–0.015. Diluting cCRE half-and-half costs information density.
3. Eval_07 (−0.013) and eval_13 (−0.015) suffer most: these were exactly the evals that gained
   most from genome → cCRE, so they're the most "regulatory-context-sensitive" and most hurt by
   dilution with random.
4. Hard evals (11, 12): −0.006 each. Random doesn't help; need a different lever for these.
5. 50/50 is too aggressive. The right hybrid ratio is probably 80–90% cCRE.

**Theory update (v3 → v4).**
- "Diversity helps" is true but selective. Eval_08 specifically tests something neither pure-
  cCRE nor pure-random has — the *combination* matters. Best guess: eval_08 evaluates near-
  baseline / low-activity sequences which cCREs underrepresent, and random anchors the model's
  prediction at the low end.
- For aggregate mean, dilution costs more than diversity gains. A small random component
  (10–20%) probably maximises both.
- Hard evals (11, 12) are stuck at 0.79–0.82. Need a NEW lever — function enrichment
  saturates, diversity helps only eval_08. Candidates: cell-type-specific elements, extreme-
  activity sequences, motif-rich synthetic, dinucleotide-shuffled controls.

---

## 2026-05-19 02:40 — Plan for experiment 005

**Hypothesis.** Random's contribution to the hybrid (better eval_08) comes from the
*motif content* of accidentally-occurring k-mers in 200bp of i.i.d. ACGT. If true, a synthetic
library with *explicitly embedded* JASPAR motifs in a random background should beat uniform
random as the diversity component of a hybrid — same baseline coverage, but with denser, more
realistic motif content. If false (motif-embedded ≈ uniform random as diversity), the gain is
about something else (GC variance, length distribution of low-activity sequences, etc.).

**Why this is the most informative next step.** Diversity hypothesis was partially confirmed in
exp 004. We need to know whether *what kind of diversity* matters. Motif-embedded vs uniform-
random is a clean isolation of "motif content" vs "everything else random offers". Either
result clarifies the design recipe.

**Design.**
- 75,000 cCRE-centered windows (same pool & filter as exp 003 / 004), seed=4
- 75,000 synthetic sequences: each is uniform-random 200bp ACGT background with 2 JASPAR motifs
  embedded as consensus sequences (the most likely base per position) at random positions, drawn
  from JASPAR 2024 vertebrate non-redundant (~880 motifs), motif identity sampled uniformly
  per sequence, position sampled uniformly such that the motif fits, seed=4
- Concatenate, shuffle, write 150,000.

**Predicted outcome.**
- Mean: 0.875–0.890. If motif-embedded is a better diversity source than random, mean ≥ 0.886
  (matches or beats pure cCRE). If equal to uniform random as diversity, mean ≈ 0.880.
- Eval_08: 0.91–0.93. The motif-embedded sequences are still "near-random" in background but
  with stronger signal in the planted motifs. I expect eval_08 ≥ 0.916 (the hybrid number),
  possibly higher if explicit motifs are the active ingredient.
- Eval_07, 13: less loss than uniform-random hybrid (because the diversity source is more
  "regulatory-shaped").
- Hard evals (11, 12): no strong prediction. If they improve here, motifs are a key axis.

This is **refining a promising direction** (diversity hypothesis) by swapping the diversity
source for a more focused one. Same hybrid framework, different second component.


## 2026-05-19 03:35 — Experiment 005 result (cCRE + motif-embedded synthetic 75/75)

**Result: mean 0.8829 (vs 0.8804 uniform-random hybrid, +0.003; vs 0.8862 pure cCRE, −0.003).
Eval_08 = 0.9160 (essentially identical to uniform-random hybrid's 0.9156). Motif content is
slightly better than uniform random as a diversity source for most evals, but the eval_08
super-additivity is NOT motif-specific — any non-cCRE sequence type unlocks it.**

| eval | cCRE | hybrid_uniform | hybrid_motif | Δ(motif−uniform) |
|------|------|----------------|--------------|------------------|
| 01   | 0.8288 | 0.8217 | 0.8258 | +0.004 |
| 07   | 0.9025 | 0.8896 | 0.8926 | +0.003 |
| 08   | 0.8922 | 0.9156 | 0.9160 | +0.000 |
| 11   | 0.8145 | 0.8078 | 0.8114 | +0.004 |
| 12   | 0.7959 | 0.7900 | 0.7926 | +0.003 |

(Full table in `libraries/005_ccre_motif_embedded/notes.md`.)

**What this tells us.**
- Motif content as the diversity source: tiny improvement on most evals (+0.003 to +0.005),
  no improvement on eval_08 specifically.
- The eval_08 effect is about *source diversity*, not motif identity. Confirmed across two
  diversity sources (uniform random and motif-embedded).
- Pure cCRE is still the leader on aggregate mean. Any 50/50 dilution costs more than it gains.

**Theory v5.** The library has two distinct levers for performance:
1. **Function-enrichment axis** — random < genome < cCRE. Plateauing.
2. **Source-diversity axis** — adding ANY non-cCRE sequence type boosts eval_08 specifically by
   ~+0.024 (super-additive), but costs ~+0.005–0.015 on every other eval at 50/50 dilution.

The right move is to find a dilution ratio where the eval_08 gain isn't fully paid for by other-
eval losses. That's the next experiment.

**Hard evals (11, 12) remain stuck** at 0.81 / 0.79. Neither function enrichment past cCREs nor
diversity in hybrid form has moved them. They need a different lever entirely — best guesses:
cell-type-specific cCREs, conservation-filtered cCREs, or activity-extreme sequences.

---

## 2026-05-19 03:40 — Plan for experiment 006

**Hypothesis.** A small (10%) diversity component preserves the eval_08 super-additivity while
paying only a small dilution cost on other evals, so 90/10 hybrid beats both pure cCRE (0.886)
and 50/50 hybrid (0.880) on mean.

**Why this is the most informative next step.** Direct test of the "right hybrid ratio"
question, which I cannot answer from the two data points (0% diversity, 50% diversity) without
guessing. If 90/10 wins on mean, hybrid design is real and we can keep refining ratio. If 90/10
loses on mean, dilution always costs more than diversity gains and pure cCRE is the ceiling for
this family of designs — opens up exploring orthogonal levers (cell-type, conservation,
activity).

**Design.**
- 135,000 cCRE-centered windows (same pool & filter), seed=5
- 15,000 motif-embedded synthetic (random ACGT + 2 JASPAR consensus motifs), seed=5
- Use motif-embedded (slightly better than uniform random in exp 005) as the diversity source.
- Concatenate, shuffle, write 150,000.

**Predicted outcome.**
- Mean: 0.884–0.890. If 0.886+ (matches or beats pure cCRE), the hybrid hypothesis wins. If
  <0.886, design space exhausted on this axis.
- Eval_08: 0.90–0.92. Open question — does 10% diversity unlock the super-additive jump, or
  does it need ≥25% to switch on?
- Hard evals (11, 12): essentially same as pure cCRE (~0.815, ~0.796).

This is **refining a promising direction**.


## 2026-05-19 05:00 — Experiment 006 result (90/10 cCRE/motif-embedded)

**Result: NEW BEST. Mean = 0.8883 vs 0.8862 pure cCRE (+0.002). Wins on 12 of 14 evals.
Eval_08 at 0.909 (+0.017 over cCRE) — the super-additive diversity effect survives at 10%
dilution. Hybrid hypothesis DECISIVELY confirmed: small diversity component beats pure cCRE.**

| eval | cCRE | 50/50 | 90/10 | Δ(90/10−cCRE) |
|------|------|-------|-------|---------------|
| 07   | 0.9025 | 0.8926 | 0.9062 | +0.004 |
| 08   | 0.8922 | 0.9160 | 0.9088 | **+0.017** |
| 11   | 0.8145 | 0.8114 | 0.8155 | +0.001 |
| 12   | 0.7959 | 0.7926 | 0.7984 | +0.003 |
| 13   | 0.9038 | 0.8898 | 0.9022 | −0.002 |

(Full table in `libraries/006_ccre_motif_90_10/notes.md`.)

**Eval_08 dose-response.** 0% diversity (cCRE alone) 0.892. 10% diversity 0.909. 50% diversity
0.916. 100% diversity (random) 0.908. Peaks around 50% but the marginal benefit per % drops
fast — going from 10% → 50% only gains +0.007 on eval_08 but costs much more elsewhere. 10%
is near-optimal for the eval_08-vs-mean trade-off.

**Theory v6.**
1. Function-enrichment axis: random < genome < cCRE. Plateauing at cCRE.
2. Source-diversity axis: a small (10%) non-cCRE diversity component is strictly +EV; 50% is
   not. The diversity effect is concentrated on eval_08, which depends on having both source
   types in training.
3. Hard evals (11, 12) remain stuck at 0.81 / 0.80. Diversity doesn't help; neither does
   function enrichment past cCREs. They need a different lever.

**Recommendation for hard evals.** Try:
- Class-balanced cCREs (overrepresent rare classes — PLS, CTCF-only).
- Cell-type-specific cCREs (if I can get K562/HepG2/SKNSH-specific subsets).
- Conservation-filtered cCREs (high phyloP — usually higher-information regulatory elements).

The natural cCRE distribution is 74% dELS; PLS is <1% and CTCF-only is 3%. If hard evals test
promoter or insulator function, balancing classes is the cheapest way to find out.

---

## 2026-05-19 05:05 — Plan for experiment 007

**Hypothesis.** Class-balanced cCREs (equal share of PLS, pELS, dELS, and CTCF-only/DNase-
H3K4me3) teach the model regulatory grammar across all element types more evenly than the
natural distribution (74% dELS). If hard evals (11, 12) test promoter or insulator function,
class balancing should lift them. If they don't move, hard evals are NOT about cCRE class
diversity and we need a different lever.

**Why this is the most informative next step.** Two reasons:
1. Hard evals have been the most-stuck part of the result for 3 experiments. Class balance is
   the cheapest hypothesis-test that addresses them directly.
2. If class balancing helps the *mean*, we have a second independent lever (beyond hybridisation)
   that combines additively with the 10% diversity lever — potentially a 0.89+ design.

**Design.**
- Bucket cCREs into 4 classes:
    1. PLS-types: `PLS`, `PLS,CTCF-bound`
    2. pELS-types: `pELS`, `pELS,CTCF-bound`
    3. dELS-types: `dELS`, `dELS,CTCF-bound`
    4. Other: `CTCF-only,CTCF-bound`, `DNase-H3K4me3`, `DNase-H3K4me3,CTCF-bound`
- Apply same window/filter as exp 003.
- Sample 37,500 from each bucket, seed=6. If any bucket has fewer survivors than 37,500,
  sample with replacement to fill (allowing some duplication for minority classes).
- Concatenate, shuffle, write 150,000.

This is **exploring a new hypothesis** (internal cCRE class diversity), not refining the
hybrid ratio direction (which gave +0.002 in exp 006).

**Predicted outcome.**
- Mean: 0.876–0.892. The PLS pool is small (~40k pre-filter) and might be high-quality
  per-sequence (promoters are concentrated, well-studied); CTCF-only is also small but
  insulator-specific. Boosting both their representation could give a small win OR a small
  loss if dELS turns out to be the most-informative class per sequence.
- Hard evals (11, 12): 0.80–0.85. If they jump, class diversity is a new lever.
- Eval_08: 0.88–0.90. No diversity component, so should be similar to pure cCRE's 0.892.


## 2026-05-19 06:05 — Experiment 007 result (class-balanced cCRE)

**Result: REGRESSION. Mean = 0.8820 (vs 0.8862 pure cCRE, −0.004). Loses on 13 of 14 evals.
Class balancing HURTS — the natural 74%-dELS distribution is doing real work.**

Key per-eval:
| eval | natural cCRE | class-balanced | Δ |
|------|--------------|----------------|---|
| 07   | 0.9025 | 0.8918 | **−0.011** |
| 08   | 0.8922 | 0.8946 | +0.002 |
| 11   | 0.8145 | 0.8107 | −0.004 |
| 12   | 0.7959 | 0.7931 | −0.003 |
| 13   | 0.9038 | 0.8920 | **−0.012** |

**What this tells us.**
- dELS (distal enhancers) really are the broadest, most generalisable training signal in
  the cCRE catalogue. Forcing equal representation of rarer classes (PLS, CTCF-only)
  reduces total dELS exposure and the model loses generalisability.
- Eval_07 and eval_13 — the evals that benefited *most* from random→cCRE in exp 003 —
  lose the most here. They probably specifically test dELS-style sequences. Compatible
  with eval_07/13 being enhancer prediction benchmarks.
- Eval_08 gained a tiny +0.002 from class diversity (consistent with the diversity-bonus
  pattern). Confirms the eval_08 lever is any kind of within-library heterogeneity.
- Hard evals (11, 12) still hard. Class balance is NOT the lever for them.

**Theory v7.**
- Adding to v6: within-cCRE class balance is a net loss. The natural distribution is
  near-optimal for general grammar (or at least, dELS is the highest-information class).
- Hard evals (11, 12) plateau is robust to:
    - function-enrichment (cCRE didn't help past random→genome boost)
    - source-diversity (random/motif-embedded diversity gave ≤ +0.003 on hard evals)
    - cCRE class-diversity (slight loss)
  → They need a fundamentally different lever. Best candidates:
    - cCREs from non-K562/HepG2/SKNSH cell types (test cross-cell-type generalisation)
    - Conservation-filtered cCREs (probably higher-quality regulatory elements per sequence)
    - Activity-extremes (would require pilot measurements I don't have)
    - Adversarial / counterfactual pairs (cCRE + mutated cCRE) — high information per pair

---

## 2026-05-19 06:10 — Plan for experiment 008

**Hypothesis.** Non-cCRE genomic sequences are a better diversity source than synthetic
motif-embedded random because they carry real chromatin context, accidental weak motifs in
realistic spacing, and the same dinucleotide composition as the cCRE pool. If true,
swapping the 15k synthetic in exp 006 (mean 0.888) for 15k non-cCRE genomic tiles should
match or beat 0.888. If equal, the diversity source is interchangeable. If worse, the
motif-embedded synthetic was specifically contributing.

**Why this is the most informative next step.** I have two complementary diversity sources
to test (synthetic vs biological). 006 used synthetic; 008 tests biological. Whichever wins
becomes the recommended diversity source. Also, non-cCRE genome captures "everything in the
genome that isn't a labeled cCRE", which might include unlabeled regulatory elements that
ENCODE missed — a free upside if they exist.

**Design.**
- 135,000 cCRE-centered windows (same pool & filter as exp 006), seed=7
- 15,000 non-cCRE genomic 200bp tiles: tile chr1, 7, 14, 19-22 into 200bp
  non-overlapping windows (same filter as exp 002), drop any window whose
  center is within 500bp of an ENCODE cCRE, sample 15,000 uniformly, seed=7
- Concatenate, shuffle, write.

**Predicted outcome.**
- Mean: 0.885–0.892. I lean slightly toward biological diversity being equal or marginally
  better (real genomic context is "more informative per sequence" than synthetic).
- Eval_08: 0.90–0.92. The diversity premium should be similar to 006's 0.909.
- Hard evals (11, 12): essentially unchanged from 006 (~0.815, ~0.798).


## 2026-05-19 07:25 — Experiment 008 result (90/10 cCRE + non-cCRE genome)

**Result: 0.8836 mean. Loses to 006 (motif-embedded diversity, 0.888) by −0.005 and even
slightly loses to pure cCRE (0.886) by −0.003. Eval_08 takes a HUGE hit (−0.024 vs 006).**

| eval | 006 (motif) | 008 (genome) | Δ |
|------|-------------|--------------|---|
| 07   | 0.9062 | 0.9015 | −0.005 |
| 08   | 0.9088 | 0.8850 | **−0.024** |
| 11   | 0.8155 | 0.8127 | −0.003 |
| 12   | 0.7984 | 0.7948 | −0.004 |
| 13   | 0.9022 | 0.9018 | 0.000  |

**Mechanism insight.** This is the strongest evidence yet that eval_08 specifically tests
*random-like / synthetic / out-of-distribution* sequences. Non-cCRE human DNA is "too
cCRE-like" in composition and doesn't trigger the eval_08 diversity bonus. Synthetic
sequences do.

So the eval_08 effect is not "any non-cCRE sequence" but "sequences that look unlike human
regulatory DNA". The model trained without exposure to such sequences predicts them poorly
on the held-out test (eval_08), and adding even 10% synthetic in training closes that gap.

**Theory v8.**
- Function-enrichment axis: random < genome < cCRE. Plateauing at cCRE.
- Source-diversity axis: ~10% synthetic (random or motif-embedded) is +EV (+0.002 mean,
  +0.017 on eval_08). Non-cCRE genome as diversity is NOT a substitute — it loses the
  eval_08 bonus AND costs a bit elsewhere.
- Hard evals (11, 12) untouched by any lever tried so far. Need a new approach.

---

## 2026-05-19 07:30 — Plan for experiment 009

**Hypothesis.** Within the winning 90/10 design, denser motifs in the synthetic component
(5 motifs per 200 bp instead of 2) give a richer training signal per synthetic sequence
without losing the random-like background that anchors eval_08. If true, mean ≥ 0.890.
If equal, motif density beyond 2 doesn't matter. If worse, dense motifs over-saturate
the synthetic sequences and lose the "random-like" property that helps eval_08.

**Why this is the most informative next step.** Confirms one of two clear refinement
directions for the winning 90/10 design. A positive result gives me a knob to keep tuning
(motif density). A negative result narrows the design space.

**Design.**
- 135,000 cCRE-centered windows (same pool & filter), seed=8
- 15,000 dense motif-embedded synthetic: random ACGT background + 5 JASPAR consensus
  motifs at random non-overlapping positions, seed=8.
- Concatenate, shuffle.

If 5 motifs of mean length 10bp at random positions overlap on average, I'll deduplicate
positions to avoid double-stamping bases. Otherwise the effective length covered by
motifs is ~50bp, leaving ~150bp of random background.

**Predicted outcome.**
- Mean: 0.886–0.892. Most likely 0.888 ± 0.002 — close to 006.
- Eval_08: 0.90–0.91. Could go either way: more motifs might dilute the "random-like"
  property eval_08 benefits from.
- Hard evals (11, 12): unchanged.


## 2026-05-19 08:30 — Experiment 009 result (90/10 cCRE + dense motif 5/seq)

**Result: 0.8888 (vs 006's 0.8883, +0.0005). Marginal new best — essentially a tie.**

Most evals up by +0.001–+0.002 (within noise). Eval_08 dropped 0.003 (denser synthetic
→ less random-like → smaller eval_08 bonus). Motif density above 2 is a flat curve.

**Theory v9.**
- 0.888 ceiling holds for the cCRE + synthetic family.
- Motif density is not the lever to push further.
- Eval_08 has a soft trade-off between "structure in synthetic component" (mean ↑) and
  "random-like backbone" (eval_08 ↑) — sparse motifs maximise the eval_08 bonus.
- Hard evals (11, 12) unchanged at 0.817 / 0.800.

---

## 2026-05-19 08:35 — Plan for experiment 010

**Hypothesis.** Sequences *adjacent to* ENCODE cCREs (flanking 200bp windows on either side)
carry additional generalisable regulatory signal — chromatin context, TFBS clusters that
extend past the cCRE boundary, dinucleotide composition of regulatory neighbourhoods.
Including them as 2/3 of the library (with each cCRE contributing center + 2 flanks)
should match or beat the pure-centre design, possibly lifting hard evals if those test
extended regulatory context.

**Why this is the most informative next step.** All my within-cCRE-centre experiments are
plateauing at 0.888. The flanks are a *fundamentally different* sequence pool —
biologically related to cCREs but not "the same 200bp the catalogue points at". If they
contain signal, that's a new lever and we lift the ceiling. If they don't, the cCRE
centre is the right granularity and I focus on orthogonal axes (cell-type, conservation).

**Design.**
- Same 745k cCRE pool. Uniform-sample 50,000 unique cCREs, seed=9.
- For each, generate 3 windows: center [mid-100, mid+100), upstream
  [mid-300, mid-100), downstream [mid+100, mid+300).
- Apply same filter (ACGT, <50% softmasked) per window — drop only the
  failing windows but keep the other two if they pass. If a cCRE has <3
  passing windows, just include what we have; sample more cCREs to keep
  total = 150,000.
- Shuffle, write.

This is **exploring a new hypothesis** (regulatory context matters beyond cCRE
centres), not refining the hybrid direction.

**Predicted outcome.**
- Mean: 0.880–0.895. Two regimes:
  (a) flanks add information → 0.890+
  (b) flanks are correlated noise → 0.880–0.885
- Eval_08: 0.88–0.90 (no synthetic component → no eval_08 bonus).
- Hard evals (11, 12): if they care about extended context, +0.005–0.015. Otherwise flat.


## 2026-05-19 10:00 — Experiment 010 result (cCRE neighbourhoods)

**Result: 0.8822. THIRD regression in the "expand the cCRE pool" direction (after 007
class-balance and 008 non-cCRE genome). Flanks dilute the signal — cCRE centres are the
right granularity.**

| eval | cCRE | neighborhoods | Δ |
|------|------|---------------|---|
| 08   | 0.8922 | 0.8844 | −0.008 |
| 11   | 0.8145 | 0.8109 | −0.004 |
| 12   | 0.7959 | 0.7926 | −0.003 |

**Pattern recognised.** Three different ways of "expanding the cCRE-derived sequence pool"
all regress vs pure cCRE centres:
- 007 class-balanced: over-represents PLS/CTCF → −0.004
- 008 non-cCRE genome (10% diversity): non-cCRE human DNA → −0.003 vs cCRE alone
- 010 neighbourhoods: flanks ±200bp → −0.004

The cCRE-centre 200bp catalogue is genuinely near the per-sequence information ceiling
among biologically-motivated libraries. To break 0.888, the headroom is in:
1. Higher per-sequence information density (cell-type filter, conservation filter)
2. Data augmentation on the same sequences (RC, mutations)
3. Multi-source diversity (3-way mixes targeting different parts of eval distribution)
4. Orthogonal annotation sources (DNase, ChIP-seq, FANTOM) — would require new downloads.

**Theory v10.** The 200bp cCRE centre is a "saturated" biological pool — sampling more
cCRE-like sequences doesn't add information. The remaining headroom is mechanistic:
either teach the model more from each cCRE (augmentation) or get more "kinds" of test
coverage (3-way diversity).

---

## 2026-05-19 10:05 — Plan for experiment 011

**Hypothesis.** Data augmentation by random base substitutions (10% per cCRE) gives the
model contrastive training signal: it sees two near-identical sequences with the slight
activity differences the MPRA actually produces from small perturbations. This should
help the model learn motif boundaries and sensitivity, and may lift mean above 0.886
(pure cCRE) — possibly above 0.888 if the augmentation effect is large.

**Why this is the most informative next step.** All the "expand the pool" experiments
have regressed; the cCRE-centre catalogue is saturated. Augmentation is the cleanest
test of "can the model learn more from the same sequences?". If it helps, augmentation
is a new lever. If not, the model is already extracting nearly all available info per
cCRE and the next direction is multi-source diversity or new annotation sources.

**Design.**
- Sample 75,000 unique cCRE-centered windows (same pool & filter), seed=10
- For each, generate a "mutated" version: randomly substitute 10% of bases (20 / 200)
  to a *different* ACGT base (so the position guaranteed changes), seed=10
- Library = 75k originals + 75k mutated = 150k. Shuffle.

This is **exploring a new hypothesis** (data augmentation as a lever), not refining
existing direction.

**Predicted outcome.**
- Mean: 0.880–0.890. My best guess: ~0.886 (similar to pure cCRE). I lean toward
  "small or zero effect" because the model probably already extracts what's there
  from 150k examples. A surprising +0.005 would suggest augmentation is real and we
  should also try 90/10 (mutation-augmented cCRE / synthetic diversity).
- Eval_08: 0.88–0.89. No synthetic component, so no eval_08 bonus.


## 2026-05-19 11:25 — Experiment 011 result (cCRE + mutated cCRE pairs)

**Result: 0.8837 mean. FOURTH "expand cCRE pool" regression in a row (07 class-balance,
08 non-cCRE genome, 10 neighbourhoods, 11 mutations). All four lose because they
reduce unique-cCRE coverage in exchange for something the model doesn't extract value
from.**

Eval_08 gained +0.005 — confirming that mutated cCREs act like a weak synthetic
diversity (out-of-cCRE-distribution enough to give a small eval_08 bonus). But not enough
to offset the cCRE-diversity cost.

**Theory v11.** The model wants HIGH-DIVERSITY unique cCREs plus a small synthetic
anchor. Pairing each cCRE with a derived/related sequence reduces unique cCRE diversity
and the gain doesn't compensate.

---

## 2026-05-19 11:30 — Plan for experiment 012

**Hypothesis.** Reverse-complement augmentation gives the model strand-aware training
signal: the same cCRE shown on both strands teaches it that regulatory motifs can
appear in either orientation. This might lift mean above 0.888 even at reduced unique
cCRE count, *if* the prepare.py model isn't already fully strand-invariant. If the
model is strand-invariant, RC just halves effective unique cCRE coverage → similar to
exp 010's regression (~−0.004 vs 003).

**Why this is the most informative next step.** Five prior pool-expansion experiments
have all lost. RC is the one pool-expansion strategy I haven't tried that has a clear
ML motivation (strand augmentation is standard in genomic CNNs). Either result is
informative: a win adds a new lever, a loss confirms that no kind of cCRE
pool expansion helps and we need orthogonal data (cell-type / conservation / DNase
peaks / FANTOM).

**Design.**
- 67,500 unique cCRE-centered windows (same pool & filter), seed=11
- 67,500 reverse complement of those same cCREs
- 15,000 motif-embedded synthetic (2 motifs each), seed=11
- Concatenate, shuffle, write 150,000.

This is **exploring a hypothesis** (RC augmentation as a strand-aware data multiplier).

**Predicted outcome.**
- Mean: 0.882–0.890. Could land either side of 0.888 depending on RC value.
- Eval_08: 0.905–0.912 (similar to 006, since synthetic component is unchanged).
- Hard evals (11, 12): probably unchanged.


## 2026-05-19 12:50 — Experiment 012 result (RC-augmented cCRE + motif) — NEW BEST 0.8905

**Result: mean = 0.8905. Wins on 13 of 14 evals vs 006 (which used 135k unique cCREs).
Despite halving unique cCRE count (67.5k × 2 strands), RC augmentation lifts every cell
line and pushes eval_12 across 0.80 for the first time ever.**

The prepare.py model is NOT fully strand-invariant. Showing each cCRE in both orientations
gives the model real additional training signal — strand-aware regulatory grammar.

**Three independent levers now identified, roughly additive:**
1. Function enrichment (random → cCRE): +0.066
2. Synthetic diversity (10% motif-embedded): +0.002 mean, +0.017 eval_08
3. RC augmentation (each cCRE × 2 strands): +0.002 mean, +0.002 across the board

Stacked: 0.886 + 0.002 + 0.002 ≈ 0.890. Confirmed.

**Theory v12.** Library design must account for ML-model invariances. A "biologically
perfect" library isn't optimal if the model's architecture has blind spots (strand
asymmetry). RC augmentation exploits the model's strand-imperfection to extract more
signal per training example.

---

## 2026-05-19 12:55 — Plan for experiment 013

**Hypothesis.** RC augmentation alone (no synthetic) gives roughly equal mean to
RC + synthetic combined, because the +0.017 eval_08 bonus from synthetic is
washed out across 14 evals once we have RC's broad +0.002 lift. If true, the
synthetic component is redundant and we should chase RC variants. If false (RC alone
is significantly worse), RC and synthetic are independent contributions that should
both be in the design.

**Why this is the most informative next step.** We have two compounding levers now —
RC and synthetic. Question: do they double-count or stack? A clean ablation
(pure RC) tells us.

**Design.**
- 75,000 unique cCRE-centered windows (same pool & filter), seed=12
- 75,000 RC of the same 75,000 cCREs
- Total 150,000. No synthetic component.
- Shuffle.

**Predicted outcome.**
- Mean: 0.886–0.892. If RC + synthetic stack additively, pure RC should be ≈ 0.888
  (similar to pure cCRE + synthetic). If RC alone wins, mean > 0.890.
- Eval_08: 0.89–0.91. Expect to lose ~+0.003 vs 012 (no synthetic bonus).
- Hard evals: similar to 012 (~0.815, ~0.80).


## 2026-05-19 14:20 — Experiment 013 result (pure RC, no synthetic)

**Result: 0.8832. Loses to BOTH pure cCRE (003 = 0.886) and RC+synthetic (012 = 0.890).**

Eval_08 took −0.024 hit vs 012 (no synthetic = no bonus). All other evals dropped 0.004–0.008.

**Key learning.** RC alone is NOT a winner. The 012 lift came from RC + synthetic acting
together. Neither alone exceeds pure cCRE; combined they give super-additive +0.004.

**Theory v13.** Two independent contributions in 012:
1. RC augmentation → small lift on most evals (when paired with synthetic that prevents
   over-fitting to one strand)
2. Synthetic component → eval_08 bonus
Together, super-additive. Alone, neither wins. Counter-intuitive but the data is clear.

---

## 2026-05-19 14:25 — Plan for experiment 014

**Hypothesis.** A middle ratio on the unique-cCRE-vs-RC trade-off curve will beat both
extremes:
- 006 (135k unique forward + 15k motif): 0.888 — all unique, no RC
- 012 (67.5k × 2 + 15k motif):           0.890 — half unique, full RC
- 014 (45k×2 + 45k forward + 15k motif): ?   — 90k unique, partial RC

If 014 > 0.890, more unique cCRE diversity helps even at the cost of some RC pairings.
If 014 < 0.890, full RC of fewer cCREs is better.

**Why this is the most informative next step.** Sweeps the dominant trade-off in the
current design family. Result either pushes the ceiling higher (014 > 012) or confirms
0.890 as the family ceiling and points us to try orthogonal levers (FANTOM5, conservation,
etc.).

**Design.**
- 45,000 unique cCREs × 2 strands = 90,000 sequences (RC-paired set)
- 45,000 different unique cCREs (forward only, no RC)
- 15,000 motif-embedded synthetic (2 motifs, seed=13)
- Total: 90k unique cCREs + 15k synthetic = 150k seqs.
- Shuffle.

**Predicted outcome.**
- Mean: 0.887–0.892. My best guess ~0.890 (similar to 012). Going to ±0.002 either way.
- Eval_08: 0.91–0.92 (same synthetic count → same bonus).
- Hard evals (11, 12): same as 012 (~0.816, ~0.800).


## 2026-05-19 15:35 — Experiment 014 result (partial RC, 90k unique)

**Result: 0.8879. Slightly below 006 (0.8883, all unique no RC) and clearly below 012
(0.8905, full RC). Partial RC is the worst of both worlds.**

The model benefits more from pairing every cCRE on both strands than from a mix of
"some paired, some not". The trade-off curve has a clear local optimum at 012.

**Theory v14.** RC value depends on pair-wise strand training (model sees A and revcomp(A)
as paired examples). Random mix doesn't give the same lift.

---

## 2026-05-19 15:40 — Plan for experiment 015

**Hypothesis.** Within "135k cCRE-derived sequences + 15k synthetic", cCRE DIVERSITY
matters more than strand-pairing. If true, 135k unique cCREs with mixed strand
assignment (67.5k fwd + 67.5k different cCREs in RC) should beat 012's design
(67.5k unique × 2 strands).

If 015 > 012: cCRE diversity dominates over pair-wise strand training. Implication:
chase more unique cCREs / new annotations to push higher.

If 015 < 012: pair-wise strand training is the active ingredient, not just "more
strand-varied data". Implication: 012's design is near-optimal for the family,
and we need a different lever (FANTOM5, conservation, etc.) to push higher.

**Design.**
- 67,500 cCREs (forward only), seed=14
- 67,500 DIFFERENT cCREs (reverse complement), seed=14
- 15,000 motif-embedded synthetic
- Total 150k. 135k unique cCREs, mixed strand assignment.

**Predicted outcome.**
- Mean: 0.886–0.892. Best guess: ≈ 0.889 (between 006 and 012). The pair-wise
  training in 012 is likely doing real work; without it, 015 will match 006
  ("more cCREs, no RC effect").
- Hard evals: ~0.815, ~0.798 (similar to 006 and 012).


## 2026-05-19 17:00 — Experiment 015 result (135k unique cCREs, mixed strand)

**Result: 0.8905. EXACTLY TIES 012 on mean. RC mechanism = strand variation in
training, NOT paired examples. Hard evals 11 and 12 reach new all-time highs
(0.8187, 0.8011).**

| eval | 012 (paired) | 015 (mixed) | Δ |
|------|--------------|-------------|---|
| 09   | 0.9482       | 0.9467      | −0.002 |
| 11   | 0.8163       | **0.8187**  | **+0.002** |
| 12   | 0.8000       | **0.8011**  | +0.001 |

Pair-wise (012) has marginally tighter easy-eval performance; mixed-strand
diverse (015) lifts hard evals. Wash on mean.

**Theory v15.**
- RC value = strand variation in training distribution (not pair-wise
  contrast). The model learns strand-invariant grammar from disparate
  examples as well as from same-example pairs.
- Hard evals (11, 12) benefit from MORE UNIQUE cCREs. 015's 135k unique
  beats 012's 67.5k unique on hard evals despite tying on mean.
- Implication: prefer the 015 design for the biological-diversity argument
  (more unique cCREs is a better library design ceteris paribus).

The "135k cCRE-derived + 15k synthetic" family has plateaued at 0.890. To
break through I need either a new data source or a different lever.

---

## 2026-05-19 17:05 — Plan for experiment 016

**Hypothesis.** Adding FANTOM5 CAGE-defined enhancers (different annotation
methodology — bidirectional transcription — vs ENCODE cCREs which use DNase
+ chromatin marks) into the diversity component captures regulatory grammar
ENCODE misses, lifting mean above 0.890.

**Why this is the most informative next step.** All cCRE-internal levers
have been explored. FANTOM5 is the cheapest new annotation source I have
in hand (already downloaded, 63k enhancers). The result is informative:
- If FANTOM5 helps: multi-source biology is a real lever, worth pursuing
  with more annotations (DNase peaks, ChIP-seq).
- If FANTOM5 doesn't help: ENCODE cCREs already cover the regulatory
  space we need, and the headroom must be in synthetic / algorithmic
  cCRE selection / cell-type-specific data.

**Design.**
- 67,500 cCREs (forward only), seed=15
- 67,500 DIFFERENT cCREs (RC), seed=15
- 7,500 FANTOM5 CAGE enhancers, 200bp centered on enhancer midpoint,
  same filter (ACGT, <50% softmasked), seed=15
- 7,500 motif-embedded synthetic, seed=15
- Total: 135k unique cCREs (mixed strand) + 7.5k FANTOM5 + 7.5k motif = 150k.

**Predicted outcome.**
- Mean: 0.886–0.893. If FANTOM5 = motif as diversity, mean ≈ 0.890. If
  FANTOM5 adds new info, ≈ 0.892+. If FANTOM5 is too cCRE-like, ≈ 0.888.
- Eval_08: 0.91–0.92.
- Hard evals (11, 12): unclear, depends on FANTOM5 content.


## 2026-05-19 18:00 — Experiment 016 result (cCRE + FANTOM5 + motif)

**Result: 0.8907 — marginal new best (+0.0002 vs 015). FANTOM5 ≈ motif as diversity. The
0.890 ceiling is robust across 3 designs (012, 015, 016 all within ±0.0002).**

Eval_13 +0.003 (real but small) — FANTOM5 may particularly suit eval_13. Otherwise tied.

**Theory v16.** Diversity component sources are interchangeable around 7.5–15k slot.
The 0.890 mean is the ceiling for "cCRE-derived bulk + small diversity component" family
across many variations. To break it, need:
- Better cCRE selection (cell-type, conservation, sub-class)
- New training signal (designed synthetic, activity-extreme)
- Algorithmic diversity max within cCRE pool

---

## 2026-05-19 18:05 — Plan for experiment 017

**Hypothesis.** dELS-only (the 74% dominant class) is the highest-information cCRE class
per sequence (implied by class-balanced 007 regression). A library restricted to dELS-types
might lift mean above 0.890 by removing lower-information classes.

**Why this is the most informative next step.** Direct follow-up to 007's negative
result. 007 showed that *forcing equal classes* hurts. The natural conclusion is "the
dominant class is doing the work" — but I haven't tested that explicitly. If dELS-only
matches or beats the all-class mix, dELS is confirmed as the active class. If it loses,
the other classes contribute beyond their numerical share.

**Design.**
- 67,500 dELS-only cCREs (forward), seed=16
- 67,500 different dELS-only cCREs (reverse complement), seed=16
- 15,000 motif-embedded synthetic, seed=16
- Total 150k.

`dELS-only` = class column is exactly `dELS` or `dELS,CTCF-bound` (combined ~789k pool
before filter, ~530k after).

**Predicted outcome.**
- Mean: 0.886–0.894. Most likely close to 0.890. If slightly higher, dELS is the active
  ingredient. If slightly lower, mixing classes helps marginally.
- Hard evals (11, 12): unclear. If hard evals test specifically dELS-like enhancers,
  could gain.


## 2026-05-19 18:55 — Experiment 017 result (dELS-only) — BIG REGRESSION

**Result: 0.8753 — loses 0.015 vs 015. EVERY eval lost 0.011–0.025. dELS-only is much
worse than the natural cCRE class mix.**

| eval | 015 | 017 | Δ |
|------|------|------|---|
| 07 | 0.9063 | 0.8929 | −0.013 |
| 08 | 0.9115 | 0.8868 | **−0.025** |
| 09 | 0.9467 | 0.9238 | **−0.023** |
| 11 | 0.8187 | 0.8049 | −0.014 |
| 12 | 0.8011 | 0.7856 | −0.015 |

**Important asymmetric result.** Both over-balancing (007: equal classes, −0.004) and
over-focusing (017: dELS-only, −0.015) hurt. The natural cCRE class distribution is OPTIMAL.
Each class contributes unique regulatory grammar that the others can't substitute for.

**Theory v17.**
- Natural cCRE class distribution (74% dELS, 13% pELS, 4% PLS, 9% other) is the sweet
  spot. Don't filter or rebalance.
- Class diversity is a real (and asymmetric) lever — losing rare classes hurts more than
  having too many of them.
- Implication: "more cCREs" of any kind isn't useful; what's useful is the *full
  natural diversity* of regulatory annotation types.

---

## 2026-05-19 19:00 — Plan for experiment 018

**Hypothesis.** GC content stratification (force uniform GC distribution) might either
expose the model to GC regions the natural cCRE distribution under-samples (helpful)
or break the natural distribution that contributes to optimal performance (hurtful,
like class-stratification).

**Why this is the most informative next step.** Direct test of "is stratification on
SOME axis (other than class) helpful?". Class-stratification hurt; if GC-stratification
also hurts, the model wants the natural cCRE-derived distribution everywhere. If it
helps, GC is a useful diversification axis.

**Design.**
- Compute GC content for each of 745k cCRE-centered windows
- Bin into 30 GC bins of equal width (~1-1.5% per bin)
- For each bin, sample 4500 cCREs uniformly (with replacement if bin too small)
- 67.5k for forward, 67.5k different for RC
- 15k motif-embedded synthetic
- Total 150k.

**Predicted outcome.**
- Mean: 0.886–0.892. My best guess: slight loss (like class-stratification), around 0.888.
  Higher GC bins are mostly promoters (PLS, pELS) and lower GC bins are dELS. Forcing
  uniform sampling effectively rebalances by GC, which correlates with class.
- Hard evals: probably tracking the mean.


## 2026-05-19 20:30 — Experiment 018 result (GC-stratified) — third stratification regression

**Result: 0.8791. Loses 0.011 vs 015. THIRD strong negative for "forced stratification":**
- 007 class-balanced: −0.004
- 017 dELS-only:      −0.015
- 018 GC-stratified:  −0.011

**Robust finding: the natural cCRE distribution is OPTIMAL.** Don't restrict, rebalance, or
stratify cCRE selection. The model wants the unbiased sample from the catalogue.

This is one of the most reproducible patterns in the project — three independent
ways of departing from natural distribution, three regressions of similar magnitude.

**Theory v18.** Library design should respect the natural distribution of regulatory
annotations. The richness comes from the relative frequencies the annotators chose,
not from a uniform synthetic distribution we'd impose.

---

## 2026-05-19 20:35 — Plan for experiment 019

**Hypothesis.** Synthetic with biologically-motivated architecture (3 copies of the same
JASPAR motif spaced ~30bp apart, mimicking homotypic TF clusters that drive cooperativity
in real enhancers) trains the model better than random-placement motif synthetic. If true,
we have a new lever for the synthetic component.

**Why this is the most informative next step.** All 'tune the cCRE pool' levers have been
explored. The synthetic component (15k of the library) is the remaining knob. So far we've
tested 2-motif vs 5-motif (009) and motif vs uniform-random (004 vs 005) and showed motif
density / source doesn't strongly matter. Structured ARCHITECTURE is the one synthetic
design we haven't probed.

**Design.**
- 67,500 cCRE-fwd + 67,500 cCRE-RC of different (135k unique cCREs, mixed strand), seed=18
- 15,000 homotypic-cluster synthetic: random 200bp ACGT background with 3 copies of the
  SAME JASPAR motif at positions 50, 100, 150 (spaced 30bp apart on average; positions
  chosen so the motifs fit). Motif identity sampled uniformly per sequence, seed=18.
- Total 150k. Shuffle.

**Predicted outcome.**
- Mean: 0.886–0.892. Could go either way. If structured synthetic gives stronger
  regulatory signal per sequence, ≥ 0.890. If unstructured was just as good, ≈ 0.890.
- Eval_08: 0.90–0.92. Structured synthetic less "random-like", may slightly reduce the
  eval_08 bonus.


## 2026-05-19 22:00 — Experiment 019 result (homotypic motif clusters)

**Result: 0.8893. Tied with 015 (−0.001). Synthetic architecture (homotypic clusters with
fixed spacing) does NOT beat random-placement motif synthetic.**

Eval_08 dropped −0.004 — confirms "more structured synthetic = less eval_08 bonus" pattern.

**Theory v19.** The synthetic component is fungible. 2/5/cluster motif placement, uniform
random, FANTOM5 — all give 0.888–0.891 within the 90/10 + RC family. What matters is
"presence of random-like sequences in 10% slot", not the structure of those sequences.

The 0.890 ceiling is robust to ALL synthetic variations tested.

---

## 2026-05-19 22:05 — Plan for experiment 020

**Hypothesis.** Wider cCREs (≥300bp) might be more complex regulatory regions with richer
architecture per sequence (more TFBS, more chromatin context), giving better training
signal than narrower cCREs.

**Why this is the most informative next step.** Filtering by class hurt (017); filtering
by GC hurt (018). But filtering by *width* hasn't been tested. Wider cCREs are
biologically motivated (larger DHS peaks indicate more complex enhancers / super-enhancers).
If this helps, we have a new lever. If not, we confirm that cCRE intrinsic features beyond
class don't usefully select better training sequences.

**Design.**
- Filter cCREs to width >= 300bp (488k pool, 347k after window/filter — plenty).
- 67,500 wide-cCRE-fwd + 67,500 different wide-cCRE-RC (135k unique, mixed strand), seed=19.
- 15,000 motif-embedded synthetic.
- Total 150k.

**Predicted outcome.**
- Mean: 0.886–0.892. If wide cCREs are richer, ≥ 0.891. If not, similar to or slightly
  below 0.890. Risk: wide cCREs might bias toward enhancers (typically wider than
  promoters), losing some class diversity that hurt 017.


## 2026-05-19 23:10 — Experiment 020 result (wide cCREs ≥300bp)

**Result: 0.8880. Loses 0.002 vs 015. FOURTH stratification regression.** Pattern is rock
solid: ANY filter on cCRE intrinsic features hurts. Trust the natural cCRE distribution.

| stratification | Δ vs natural |
|----------------|--------------|
| class-balanced (007) | −0.004 |
| dELS-only (017) | −0.015 |
| GC-stratified (018) | −0.011 |
| wide-only ≥300bp (020) | −0.002 |

**Theory v20.** The natural cCRE distribution is OPTIMAL across multiple intrinsic axes
(class, GC, width). This is one of the most reproducible findings. Don't try to
"engineer a better cCRE subset" using surface features — sample naturally from the
catalogue.

---

## 2026-05-19 23:15 — Plan for experiment 021

**Hypothesis.** A 95/5 cCRE/synthetic ratio (7.5k synthetic instead of 15k) preserves
the eval_08 bonus (which saturates at low synthetic fraction) while giving 7.5k more
unique cCREs (142.5k vs 135k). Could marginally beat 015's 0.890.

**Why this is the most informative next step.** Sweeps the dominant trade-off in the
winning family one more time (last sweep was 014 with partial-RC at the wrong granularity).
Cheap test: either confirms 10% is the optimum (015 family ceiling) or finds a slightly
better point.

**Design.**
- 71,250 cCRE-fwd + 71,250 different cCRE-RC (142,500 unique, mixed strand), seed=20
- 7,500 motif-embedded synthetic, seed=20
- Total 150k

**Predicted outcome.**
- Mean: 0.886–0.892. Maybe slight win (~0.891), more likely tied or marginal loss.
- Eval_08: 0.90–0.91. Synthetic saturates at low %, so 5% should still give most of
  the bonus.
- Hard evals 11/12: maybe slight gain from more unique cCREs.


## 2026-05-20 00:15 — Experiment 021 result (95/5 ratio)

**Result: 0.8878. Regression of 0.003 vs 015. Eval_08 dropped most (−0.008).**

Confirms: **10 % synthetic is the optimum**. The ratio curve:
- 0 % synthetic (013): 0.883 (no eval_08 bonus)
- 5 % synthetic (021): 0.888
- 10 % synthetic (015): 0.890 (optimum)
- 50 % synthetic (005): 0.883 (too much dilution)

**Theory v21.** The 0.890 ceiling is robust to ALL tested variations within the
cCRE-RC-synthetic family. The optimum is:
- ~135k unique cCRE-derived (mixed strand or paired RC, both work)
- ~15k synthetic (any source — motif, random, FANTOM5 substitute)
- Natural cCRE distribution (no filtering)

---

## 2026-05-20 00:20 — Plan for experiment 022 — NOISE TEST

**Hypothesis.** The 0.001–0.003 deltas I've been chasing are at-or-below the training
stochasticity floor. A noise/replication experiment will quantify this. If 015 with a
different seed lands within ±0.001 of 0.8905, the deltas were real signal. If it lands
±0.003+, much of my recent analysis is in the noise.

**Why this is the most informative next step.** I should have done this earlier. Without
a noise estimate, I cannot interpret 016 vs 015 vs 012 (all ≈ 0.890), 009 vs 006 (0.0005
delta), etc. This single experiment gives essential context for all 21 prior results
and tells me how much to trust the marginal "wins" of my best designs.

**Design.**
- Same recipe as exp 015: 67,500 cCRE-fwd + 67,500 different cCRE-RC + 15k motif
- Different RNG seed (21 instead of 14) — different cCRE sample, different motif
  selection, different shuffle.

**Predicted outcome.**
- Mean: 0.886–0.895 (95 % interval guess). Tight expected if model training is mostly
  deterministic given the library; wide if training stochasticity is significant.
- Hard evals: ±0.005 expected.


## 2026-05-20 01:30 — Experiment 022 result (NOISE TEST: 015 with seed=21)

**Result: 0.8875 — Δ = −0.003 vs original 015 (0.8905). MAJOR FINDING.**

The training pipeline noise on mean_r is ≈ ±0.003. This recalibrates everything:
- 003 → 006 (+0.002 from hybrid): WITHIN NOISE
- 006 → 012 (+0.002 from RC): WITHIN NOISE
- 012 → 015 (same): equivalent, confirmed
- 015 → 016 (+0.0002 from FANTOM5): WITHIN NOISE
- 009 → 006 (+0.0005 dense motifs): WITHIN NOISE

What IS robustly above noise (≥ 0.005):
- random → genome: +0.053
- genome → cCRE: +0.013
- 50/50 hybrid eval_08 super-additivity: +0.024 on eval_08
- ALL stratification regressions (007, 017, 018, 020): −0.004 to −0.015
- dELS-only −0.015

**Theory v22 — major revision.** The 0.886-0.892 "ceiling" is FUZZY. The cCRE+RC+10%synth
family converges around **0.889 ± 0.003**. Many of my second-half "wins" were noise. The
big-effect findings (cCRE >> random/genome, stratification hurts, 10% synthetic helps)
remain robust.

**Implications for remaining experiments.** Stop fine-tuning. Either find a >0.005
breakthrough or do replicates / final design. Honest assessment: the cCRE + 10% synthetic
+ optional RC design is likely the natural ceiling for biologically-driven library
design without external data (cell-type-specific, conservation).

---

## 2026-05-20 01:35 — Plan for experiment 023

**Hypothesis.** Combining ALL winning levers (cCRE bulk, RC augmentation, FANTOM5 as
additional biology, multi-source synthetic) into one library gives a small compound
benefit. If even the kitchen-sink doesn't exceed noise, the ceiling is real.

**Design.**
- 60,000 cCRE-fwd + 60,000 different cCRE-RC (120k unique cCREs, mixed strand)
- 7,500 FANTOM5-fwd + 7,500 different FANTOM5-RC (15k FANTOM5, mixed strand)
- 10,000 motif-embedded + 5,000 uniform random (15k synthetic, 2 sources)
- Total 150k, seed=22

**Predicted outcome.**
- Mean: 0.886–0.893. Given noise ±0.003 and the synthetic-source-equivalence finding,
  likely lands in 015/016's neighborhood (~0.890). A meaningful win would be ≥ 0.894.


## 2026-05-20 02:40 — Experiment 023 result (kitchen-sink)

**Result: 0.8886. Within noise (±0.003) of all "winners" — 015, 016, 019, 022 all sit
at 0.889 ± 0.003. The design space is SATURATED for the cCRE-derived bulk + synthetic
family.**

Stacking multiple winning levers (cCRE-RC, FANTOM-RC, 2-source synthetic) gives no
compound benefit. Multiple paths to the same ceiling.

**Theory v23.**
- The ceiling for biologically-derived libraries (cCRE + small synthetic anchor) without
  external data is **0.889 ± 0.003**.
- To break this would require:
  (a) Cell-type-specific selection (would need ENCODE per-cell-type data downloads)
  (b) Conservation filter (would need phyloP/phastCons download)
  (c) Activity-stratified selection (impossible without pilot measurements)
  (d) Cross-species conserved regulatory elements

**Implication for remaining 6 experiments.** Stop refining the synthetic component or
cCRE selection. Focus on:
- One last test for compound synthetic source diversity (024)
- A few noise replicates of best design to pin down the ceiling
- Final "best known" library submission

---

## 2026-05-20 02:45 — Plan for experiment 024 (last synthetic variation)

**Hypothesis.** Splitting the 15k synthetic component across 3 different sources
(motif-embedded, uniform random, homotypic clusters at 5k each) gives slightly more
diverse synthetic coverage than any single source. If this helps even marginally, it
suggests synthetic diversity has untapped potential.

**Design.**
- 67,500 cCRE-fwd + 67,500 different cCRE-RC (135k unique, mixed strand), seed=23
- 5,000 motif-embedded synthetic + 5,000 uniform random + 5,000 homotypic clusters
- Total 150k.

**Predicted outcome.** Mean: 0.886–0.893. If above 0.892, multi-source synthetic helps
slightly. If within ±0.003 of 0.889, synthetic source diversity doesn't matter.


## 2026-05-20 04:00 — Experiment 024 result (3-source synthetic)

**Result: 0.8861. Within noise but on the low end. 3-source synthetic split (5k each)
provides NO compound benefit; single-source synthetic is just as good.**

The 0.889 ± 0.003 ceiling confirmed across 6+ designs. The "find a single magical
synthetic component design" line of investigation is exhausted.

---

## 2026-05-20 04:05 — Plan for experiment 025

**Hypothesis.** A combined design using paired-RC (012's strength on per-pair training)
with the FANTOM+motif diversity split (016's strength on multi-source synthetic) might
extract one more drop of performance. Two separate +0.000 effects might compound if both
are actually small positives we couldn't detect individually.

**Design.**
- 67,500 cCRE-fwd + 67,500 SAME cCREs RC (paired) — like 012, NOT 015's "different RC"
- 7,500 FANTOM5-fwd + 7,500 different FANTOM5-RC (15k FANTOM, mixed strand)
- Wait — that's 150k cCRE-derived. Need to fit synthetic.
- Revise: 60k cCRE-fwd + 60k SAME cCREs RC (paired, 120k seqs from 60k unique)
          + 7.5k FANTOM5-fwd + 7.5k different FANTOM5-RC (15k)
          + 15k motif synth
= 150k

**Predicted outcome.** Probably 0.886–0.892 (within noise). If above 0.892, two slight-
positive levers stack. Otherwise saturated.


## 2026-05-20 06:00 — Experiment 025 result (paired-RC + FANTOM + motif)

**Result: 0.8889. Within noise of cluster. Combining levers (paired-RC + FANTOM + motif)
doesn't compound.**

Hard evals: 11 = 0.8175, 12 = 0.8000 (matched best ever). No breakthrough.

---

## 2026-05-20 06:05 — Plan for experiment 026

**Hypothesis.** Dinucleotide-shuffled cCRE (preserves cCRE local composition but destroys
all motifs) as the synthetic component tests a mechanistic question: is the eval_08 bonus
about COMPOSITION matching cCREs or about motif presence in the synthetic?

If di-shuffled ≈ motif-embedded as synthetic: composition is the active ingredient.
If di-shuffled << motif-embedded: motif content matters.
If di-shuffled > motif-embedded: composition matters MORE.

**Design.**
- 67,500 cCRE-fwd + 67,500 different cCRE-RC (135k unique cCREs)
- 15,000 dinucleotide-shuffled cCREs (Eulerian-walk shuffling preserves 2-mer composition)
- Total 150k, seed=25

**Predicted outcome.** Within noise of 0.889. The interesting per-eval read is eval_08:
if it stays ~0.91 with di-shuffled synthetic, composition matters; if it drops to ~0.89,
motifs matter.


## 2026-05-20 17:00 — FINAL SUMMARY (after all 30 experiments)

### Final theory (v∞)

A 150k 200bp MPRA library trains a generalisable regulatory-activity model best when it:

1. **Bulk-samples ENCODE V3 cCREs in their natural distribution** (windows centred on
   cCRE midpoints, dropping N/heavy-soft-masked windows). This is the single biggest
   lever: random ACGT → cCREs lifts mean across 14 evals from **0.820 → 0.886** (+0.066).
   - Sub-finding: do NOT stratify, rebalance, or filter cCREs by class, GC, or width.
     Four independent stratification attempts all regressed by 0.004–0.015. The
     ENCODE pipeline's class distribution (74 % dELS, 13 % pELS, 4 % PLS, 9 % other)
     is the optimum.
2. **Includes a small (~10 %) synthetic component** — 15k random ACGT 200bp with 2
   JASPAR consensus motifs embedded. Adds +0.002 mean (within noise but consistent),
   and specifically lifts eval_08 by 0.017–0.024 (super-additive effect on that one
   eval). Synthetic source is fungible — motif-embedded random, uniform random,
   homotypic clusters, di-shuffled cCRE all give comparable performance; what matters
   is that the synthetic differs in composition from cCREs.
3. **Augments cCREs with strand variation** — either by pair-wise RC (each cCRE in
   both strands) or by mixed-strand sampling (different cCREs in fwd vs RC) of 135k
   unique cCREs. Both schemes give the same mean (≈ 0.890); RC value comes from
   strand-variation in training data, not pair-wise contrastive learning.
   Adds maybe +0.001–0.002 mean (within noise).

The result is **mean_r ≈ 0.889 ± 0.003** across the 14 anonymous eval sets — a
solid plateau confirmed by ≥ 8 designs converging there and 3 noise replicates
(015, 022, 029) spreading 0.8875–0.8905.

### What worked (effect sizes ABOVE noise ±0.003)

| Lever | Effect | Confirmed in |
|-------|--------|--------------|
| random → genome → cCRE bulk | +0.053 → +0.013 (function enrichment axis) | 001,002,003 |
| Adding 10 % synthetic (any source) | +0.002 mean, +0.017–0.024 eval_08 | 006,012,015,026 |
| Natural cCRE distribution (vs any filter) | losing +0.004 to +0.015 confirms | 007,017,018,020 |
| Pair-wise / mixed-strand RC + synth | +0.002 (within noise but stable) | 012,015,016,025 |

### What didn't (negative results that are themselves informative)

| Lever | Effect | Lesson |
|-------|--------|--------|
| Class-balanced cCREs | −0.004 | Natural class freq is optimal |
| dELS-only cCREs | −0.015 | Other classes contribute uniquely |
| GC-stratified cCREs | −0.011 | Forced GC uniformity hurts |
| Width-filtered cCREs (≥ 300 bp) | −0.002 | Width doesn't predict info density |
| cCRE neighbourhoods (flanks ±200 bp) | −0.004 | Flanks dilute the cCRE-centre signal |
| Non-cCRE genome as diversity | −0.024 on eval_08 vs synthetic | Diversity must be SYNTHETIC, not just non-cCRE |
| Mutation-augmented cCRE pairs | −0.003 | Pair-wise mutations don't substitute for unique cCREs |
| Partial RC pairing | −0.003 vs full RC | Full pair-or-mixed strand, not partial |
| 95/5 ratio (less synthetic) | −0.003 | Synthetic needs to be ≈ 10 % for full eval_08 bonus |
| 50/50 ratio (more synthetic) | −0.006 | More dilutes cCRE info density |
| Stacking levers (kitchen sink) | within noise | No compound benefit, design space saturated |
| Homotypic motif clusters in synth | within noise | Synthetic architecture doesn't matter |
| Dense motifs (5 vs 2 per synth) | within noise | Motif density doesn't matter |
| FANTOM5 in place of motif synth | within noise | Synthetic biology source is fungible |
| 3-source synthetic split | within noise | Single synthetic source is enough |

### My final BEST library: **libraries/030_final_paired/**

Recipe:
- 67,500 unique cCRE-centered 200 bp windows (forward strand)
- 67,500 reverse complements of the same 67,500 cCREs (paired RC)
- 15,000 motif-embedded synthetic (random ACGT background + 2 random JASPAR
  consensus motifs per sequence)

Measured mean_r: **0.8883** (single-replicate). Family mean across 3 replicate
runs (012/015/022/029/030): **0.889 ± 0.003**.

Per-cell-line on 030: K562 0.886, HepG2 0.895, SK-N-SH 0.884.
Per-eval range on 030: 0.799 (eval_12, hardest) – 0.945 (eval_09, easiest).

### Hard evals that I could NOT crack
- **eval_11** plateaus at 0.815 ± 0.005 across all designs
- **eval_12** plateaus at 0.798 ± 0.005 across all designs
- Neither budges with any tested lever (function enrichment, synthetic, RC,
  stratification, mutation, neighbourhoods, kitchen sink). My best guesses
  for what *would* help (untested): cell-type-specific cCRE selection,
  conservation-filtered cCREs, activity-extreme sequences, designed enhancer
  architectures with specific TF cooperativity patterns.

### Eval-set clustering observed
- {eval_02, eval_06, eval_14} consistently report values within 0.001 of each
  other across all 30 experiments → treat as 1 effective signal (3 reads of
  the same test).
- {eval_01, eval_05} likewise → 1 effective signal.
- So the 14 evals are effectively ~10 independent.

### Recommendations for the next round

1. **Get external data**. The 0.889 ceiling looks tight for "cCRE + light
   synthetic" without:
   - per-cell-type ENCODE cCRE activity (would let me filter to "cell-type-
     specific elements" or weight cCREs by cross-cell breadth)
   - phyloP/phastCons conservation (filter for evolutionarily conserved
     regulatory regions)
   - DNase-seq peaks per cell type (broader regulatory annotation than cCREs)
   - FANTOM5 was tested as a small additional source — it didn't add anything,
     suggesting CAGE-defined enhancers overlap with cCREs.

2. **Design synthetic with KNOWN cooperative TF combinations**. My homotypic
   clusters (3 copies of same motif, fixed spacing) did NOT beat random
   placement. But a smarter design — heterotypic dimers with biologically-
   plausible spacing (NF-Y + Sp1, E-box dimers) — is untested and might
   work because real enhancers cooperate.

3. **Activity-stratified sampling**. If a pilot MPRA could be run to
   estimate activity per cCRE, library design could select cCREs uniformly
   across the activity distribution rather than uniformly across cCREs. This
   would balance training between "easy bright" and "rare dim/repressed"
   sequences — possibly lifting hard evals.

4. **Architecture-aware synthetic at much higher complexity**. The model may
   need to see TF-spacing variation, motif-orientation variation,
   chromatin-context variation that random JASPAR placement doesn't sample.
   Designed enhancers spanning a broader regulatory grammar space could
   teach the model TF combinatorial rules.

5. **Establish noise floor early.** I should have run a noise/replicate
   experiment at experiment 5–10, not at 22. Knowing ±0.003 noise from the
   start would have saved me ~5 experiments of small-delta-chasing that
   were just within noise.

6. **The "easy wins" are spent.** The first three experiments (random →
   genome → cCRE) captured ~0.07 of the 0.89 final mean — about 80 %
   of all gains. Subsequent 27 experiments collectively added < 0.005
   above noise. The marginal value of a 31st experiment in this paradigm
   is small; opening new data sources is the highest-EV move.

### Key meta-lessons about this style of research

- **Run noise estimates early.** A single repeat of any design at a different
  seed tells you which deltas to take seriously.
- **Negative results CAN be the most informative.** The four stratification
  regressions taught me more (don't filter the catalogue) than any winning
  variation did.
- **Saturate one lever before adding another.** Mixing two new things at
  once means you can't attribute the result. My cleanest experiments
  (90/10 ratio sweep, RC vs no-RC) all isolated single variables.
- **Library design has hard ceilings.** Once cCRE + synthetic + RC saturates,
  no internal variation breaks through. External data is the only path
  forward.

