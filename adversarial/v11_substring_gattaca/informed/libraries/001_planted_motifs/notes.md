# E1: planted_motifs

## Design
- 50k 200bp sequences. 80% with 1-4 planted canonical TF motifs (~30 TF
  families), 20% pure random 50% GC backbone.
- Single seed (0).

## Result
- mean_r averaged across 14 evals ≈ **0.808**
- eval_01: 0.8238

## Comparison vs baselines (single-seed E1 vs 5-seed baselines)
| eval | E1 | random_uniform | gc_50 | delta vs random |
|------|----|----------------|-------|---------|
| 01   | 0.8238 | 0.8566 | 0.8591 | -0.033 |
| 04   | 0.8217 | 0.8631 | 0.8667 | -0.041 |
| 07   | 0.7569 | 0.7920 | 0.7939 | -0.035 |
| 08   | 0.7656 | 0.7726 | 0.7729 | -0.007 |
| 13   | 0.7898 | 0.8213 | 0.8220 | -0.032 |

E1 LOSES on every eval by ~0.03 vs the random/gc_50 baselines.

## Interpretation
Planting explicit motifs into random sequences degrades model performance
relative to pure random. The most likely explanation: the eval set is
itself drawn from a distribution close to random 50% GC (or natural
sequences with similar k-mer statistics). Planting motifs creates a
training-eval distribution shift — the model is trained on
motif-enriched sequences but evaluated on sequences without that
enrichment. Even though motifs are "biological signal", they are
off-distribution and harm generalization to this eval set.

## Update to theory
- The eval distribution closely resembles random 50% GC. Adding explicit
  structured motifs creates covariate shift that hurts.
- "Informational value" here is dominated by matching the eval input
  distribution.
- To improve over random_uniform, need to either (a) precisely match
  higher-order statistics without injecting structure, or (b) provide
  diversity in *activity-relevant* features that exist in random sequences
  by chance.

## Next
E2 will replicate pure random 50% GC to (1) calibrate my pipeline to the
published baseline of 0.857 and (2) measure single-seed run variance so I
know how big a delta is meaningful.
