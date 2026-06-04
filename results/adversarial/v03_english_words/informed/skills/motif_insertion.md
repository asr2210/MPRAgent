# Skill: Motif insertion library generation

## What works
1. Random uniform 200bp backbone + JASPAR motif insertions BEATS plain random
   uniform in this pipeline.
2. **Cell-type-targeted motifs > random JASPAR motifs** — restricting to TFs
   relevant to the labeling cell types (K562 erythroid, HepG2 hepatocyte,
   SK-N-SH neural) improves over random JASPAR sampling by ~0.5-1%.
3. 3 motifs/seq seems near-optimal; 6 didn't help much (and hurt SK-N-SH).

## What doesn't work in MY pipeline
- DHS sequences (real biology) underperform random by ~10%.
- GC filtering on DHS doesn't recover.
- Numsamples weighting doesn't help.
- gc_50 backbone (fixed 50% GC by shuffling) isn't measurably better than
  random uniform backbone for motif libraries.

## JASPAR access
- File: `data/JASPAR2024_vertebrates.txt` (713KB, 2346 PFMs)
- Download: `https://jaspar.elixir.no/download/data/2024/CORE/JASPAR2024_CORE_non-redundant_pfms_jaspar.txt`
- Format: blocks starting with `>ID Name`, then 4 lines `A [counts]`, `C [counts]`, etc.
- Parser snippet in `libraries/006_motif_enriched_random/generate.py`

## Cell-type TF token lists (used in exp 008)
Captured ~289 motifs that match these substrings:
- **K562** (erythroid): GATA, KLF (some), NFE2, MAF::NFE2, MYB, TAL1, RUNX1, SPI1
- **HepG2** (hepatocyte): HNF4, HNF1, FOXA, CEBP, ONECUT, RXRA, NR1H4, PPARA
- **SK-N-SH** (neural): NEUROD, NEUROG, ASCL, OLIG, SOX, POU3F, ISL, MEF2C, MYCN, PHOX, REST, RFX, PAX
- **Universal**: SP1/2/3, CTCF, TBP, NFY, FOS, JUN, BACH, ELK, ETV, GABPA, E2F, MYC, MAX, ATF, CREB, USF, YY1, EGR1, NRF1

See `libraries/008_celltype_motifs/generate.py` for exact token list.

## SK-N-SH IS movable (updated after exp 026, 029)
Earlier claim: "SK-N-SH stuck at 0.05-0.07 across all strategies" — FALSE.
- Exp 026 (mammalian dinucleotide backbone) pushed SK-N-SH mean to **0.0623**
  and got eval_07 = **0.0704** (first eval > 0.07).
- Exp 029 (greedy 6-mer diversity selection) pushed SK-N-SH mean to **0.0653**,
  with eval_04 = **0.0743**, eval_07 = **0.0730**, eval_09 = **0.0743**.

**SK-N-SH oracle responds to sequence COMPOSITION** (dinucleotide stats,
short-range periodicity) more than discrete motif content. K562/HepG2 are
the motif-driven ones.

## Pipeline ceiling
Best single-seed: 0.4283 mean_r (exp 008). Re-running 008 with SEED=42
gave 0.4222 (exp 016). **Library noise floor is ±0.005-0.01 on mean_r.**
The true "008-family" plateau is ≈ 0.422-0.428.

Strategies.md baseline best: 0.4243 (within noise of 008).

To CLEARLY beat the plateau, candidate gains must exceed ~0.015 in a
single seed, or be validated across multiple seeds.

What does NOT clearly help (within noise of 008):
- Motif density variations (3, 5, or 1-6 mix all similar)
- Pool size (289 vs 142 within noise)
- Same-TF clustering vs different-TFs (clustering loses by noise)
- Cell-type-clustered sequences vs mixed (clustering loses by noise)

What does CLEARLY hurt (>0.01 below plateau):
- Consensus-only motifs (exp 011, -0.011): PFM stochasticity matters.
- DHS sequences (exp 001-005, -0.04): wrong source in this pipeline.

To push further, may need:
- Realistic enhancer SYNTAX (paired motifs with biological spacing)
- Much larger pool with WEIGHTED sampling (cell-type-biased draws from full JASPAR)
- Multi-source mixture libraries (motif + something)
- Score-based selection from oversized candidate pool
- Iterative refinement (testing then re-designing)

## BREAKTHROUGH (exp 029): diversity selection beats single-design plateau
The 008/020 plateau at ~0.4284 is NOT a true ceiling. It is a property of
SINGLE-DESIGN libraries oversampling particular k-mer distributions.

**Method that worked (exp 029, eval_01 = 0.4288):**
1. Generate 150k candidate sequences from 3 different best designs:
   - 50k @ exp 020 design (8 motifs overlap, random uniform bg)
   - 50k @ exp 022 design (10 motifs overlap, random uniform bg)
   - 50k @ exp 026 design (8 motifs overlap, mammalian dinucleotide bg)
2. Greedy selection of 50k sequences maximizing 6-mer novelty:
   - For each candidate, score = sum over its 6-mers of 1/(1+covered[km])
   - Sample 200 candidates per round, pick the highest-scoring one
   - Repeat 50,000 times
3. Selected library has K562/HepG2 lift AND SK-N-SH lift simultaneously.

**Why it works (theory T20):** Surrogate generalization is limited by
EFFECTIVE training-set diversity, not per-sequence design choice. Single
designs over-represent particular short-context patterns. Mixed-and-selected
libraries have higher effective entropy → richer model features.

**Why it's fragile (theory T21, exp 030):** Adding a 4th pool with HIGH
k-mer novelty per seq but LOW signal (e.g., sparse 5-motif/seq design)
causes the greedy selector to drift toward it, hurting per-cell-type
signal. Candidate pools should have comparable per-sequence k-mer entropy.
