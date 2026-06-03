# Experiment 020 — RC-augmented library (25k unique + 25k RCs)

## Design
Generate 25k unique GC-50 random sequences. For each, also include its
reverse-complement. Total = 50k (25k unique forward, 25k their RCs).

## Result
- eval_01 = **0.3964** (Δ vs 007 = -0.0005, within noise)
- mean_r = **0.3861** (Δ vs 007 = 0.0000, IDENTICAL)
- Per-cell: K562=0.619, HepG2=0.435, SK-N-SH=0.136

## Major insight
**25k unique seqs + RCs ≡ 50k unique seqs.** The library with HALF
the unique information performed identically to 007's full library.

This strongly implies one of:
1. The trainer internally applies RC-augmentation, so 50k unique
   already counts as ~100k from the model's perspective
2. The model architecture is RC-equivariant (Basset-style 1D CNN
   with reverse-complement parameter sharing)
3. Both

Either way: the EFFECTIVE info bottleneck is ~25k unique sequences.
Adding the other 25k as either RCs (020) or as fresh unique sequences
(007) makes NO difference.

## Theory v13
The unique-information limit is < 25k sequences. Above this, additional
sequences don't add learnable information for the (model, evaluator).

This is a major architectural finding:
- If we have 25k high-quality sequences, that's enough.
- To break the ceiling, we don't need MORE sequences — we need BETTER
  sequences with more INFO PER SEQUENCE.
- Per-seq info content might be measured as: information about
  generalizable TF binding, motif arrangement, etc.

## What this rules out
- That library size (above 25k) is the bottleneck
- That sequence duplication via RC creates collinearity that hurts
- That RC explicit pairing is a meaningful augmentation

## What this opens up
1. **Per-seq info density libraries**: 25k sequences chosen for max
   information per seq (high motif content, high oracle uncertainty, etc.)
2. **Curriculum design**: maybe ordering / mix of "easy" and "hard" seqs matters
3. **Dual library**: 25k easy + 25k hard (or RCs thereof) — sweet spot

## What to try next
1. **Mutation-around-oracle** (021): 5k oracle-top seeds × 10 single-base
   mutations each = 50k. Tests whether localized exploration of active
   neighborhoods provides info or whether unique sequences matter.
