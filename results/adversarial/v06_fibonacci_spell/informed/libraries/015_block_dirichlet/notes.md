# Experiment 015 — Block-Dirichlet (4 × 50bp blocks)

## Result
eval_01 = **0.1348** vs Dirichlet(0.5) 0.138 ± 0.003 (002+014 noise floor)

Block-Dirichlet is meaningfully worse (~-0.003). Within-sequence
discrete block variation hurts more than the smooth gradient
(009 = 0.1369).

## Pattern across within-sequence structure experiments
| structure type | result | gap vs noise floor |
|---|---|---|
| Markov-Dirichlet (003) | 0.1346 | -0.003 |
| Block 4x50bp (015) | 0.1348 | -0.003 |
| Composition-iid (008) | 0.1354 | -0.003 |
| Gradient smooth (009) | 0.1369 | within noise |
| JASPAR motifs (006) | 0.1364 | within noise |

## Lesson
Smooth structure (gradient) ≈ noise. Discrete structure (blocks,
Markov runs, real DHS) costs ~0.003 below noise floor. The evaluator
seems to prefer:
- iid composition draws over discrete structural changes
- Or: structure changes that don't break the position-independence
  assumption it learned

## Definitive: within-seq structure axis is closed
Pure iid Dirichlet(0.5) is the local maximum. No within-seq variation
helps. Going forward: seed scan + final refinement.
