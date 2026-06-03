# 003 — JASPAR motif-embedded synthetic

## What I built
50,000 sequences = 200bp random ACGT backgrounds with 4-8 JASPAR 2024 CORE vertebrate (879 motifs) embedded at random non-overlapping positions. Each motif sampled per-column from the PWM (not consensus), random strand.

## Predicted
If motifs are the core regulatory grammar (H1), eval_01 should beat dirichlet 0.0768. If under-fit at 50k (H7), eval_01 ~0.07-0.08.

## Result
- eval_01 = 0.0671 — WORSE than random_uniform (0.0711) by 0.004 and worse than DHS-stratified (0.0739) by 0.007.
- eval_07 = 0.1328 (vs 0.1454 in exp 001) — biology-sensitive eval also dropped.
- eval_08 = 0.0602 (vs 0.0627 in exp 001) — hard eval dropped.
- Worst result so far on eval_01.

## Interpretation
Embedding motifs into random backgrounds HURT performance vs both random and DHS. Three candidate explanations:
1. **Composition homogenization**: motifs themselves have specific compositions (often GC-rich, palindromes, etc.). Embedding 4-8 motifs per 200bp partially homogenizes the per-sequence composition distribution. Library compositional variance shrinks vs pure random → less learnable variance → lower fit.
2. **Spurious anti-natural structure**: random backgrounds + random motif placement create sequences that aren't on the natural data manifold. The MPRA activities are noisy and uncorrelated with the embedded motifs in a coherent way that the model can fit.
3. **Lost length flexibility**: random backgrounds are positionally uniform; embedding fixed motifs creates positional discontinuities the model misfits.

Most likely (1): compositional variance is the dominant learnable signal in this harness, and embedding structured elements REDUCES it.

## Theory updates
- H1 (motifs are core grammar): NOT supported at 50k. Motif content alone doesn't translate to learnable signal.
- H5 (dynamic range matters): supported. Embedded motifs narrowed range → score dropped.
- H7 (model under-fit at 50k, learns composition): STRONGLY supported. Composition variance is the dominant lever, not biology.

This is a coherent story emerging: in this 50k harness, the trained model can only learn coarse features (composition trends, perhaps a few CpG/GC-rich vs AT-rich associations). Motif grammar requires more data than 50k provides.

## What to try next
- EXP 004: explicit composition-variance test — Dirichlet(0.3) vs Dirichlet(1.0). If higher per-sequence composition variance improves eval_01 further, composition is confirmed as the lever.
- Then: combine composition diversity with strong-signal biology. E.g., very-high-signal DHS that are compositionally extreme.
