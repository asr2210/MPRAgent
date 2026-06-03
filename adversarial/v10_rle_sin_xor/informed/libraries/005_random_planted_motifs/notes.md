# Experiment 005 — Random with planted JASPAR motifs

## Hypothesis
Random sequences fail to reach HepG2 r > 0.56 because they lack TF binding motifs.
Planting 8 random JASPAR motifs per 200bp sequence (from 873 vertebrate TFs) should
provide rich motif content the oracle can learn from.

## Result
eval_01 = 0.5161. K562=0.984, HepG2=0.567, SK-N-SH=0.
**Slightly WORSE than synth_random (0.521).** K562 dropped 0.99 → 0.98 (motif content
slightly disrupts the K562 oracle's preferred input).

## Verdict
Motif planting from broad TF set does NOT improve HepG2 prediction. The HepG2 oracle
is not motif-driven — at least not by typical JASPAR TFBS planted in random backbones.

Combined with prior 4 experiments, the picture is clear:
- K562 oracle saturates at ~0.99 for any high-entropy random library
- HepG2 oracle saturates at ~0.56 for any 50% GC sequence content
- SK-N-SH r = 0 regardless of library

The strategies.md baselines suggest random_uniform is the strongest simple strategy
in this pipeline. Anything deviating hurts.
