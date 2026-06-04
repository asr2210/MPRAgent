# 012 — Block composition (2x100bp halves)

## What I built
50,000 sequences. Each = two 100bp blocks. Each block has its own independent
Dirichlet(0.3) composition draw.

## Result
- eval_01 = 0.0748 (vs exp 004 single-block Dirichlet(0.3) = 0.0786). **-0.004 loss**.
- eval_08 dropped 0.0716 → 0.0618 (more severe loss on noisiest eval).
- mean ≈ 0.0907 (vs 0.0954 for exp 004).

## Interpretation
Block composition HURT. The likely mechanism: averaging two independent Dirichlet(0.3)
draws pulls the per-sequence OVERALL composition closer to the middle. The model sees
less extreme overall compositions → less per-sequence composition spread → less learnable
signal.

This is consistent with the theory: **the model learns from overall per-sequence
composition, not from positional sub-structure within sequences.**

eval_08 confirmed as the most composition-spread-sensitive eval (drops sharply whenever
overall composition is moderated).

## Hypothesis killed
H12: positional block composition adds learnable info → **REJECTED**

## What to try next
Test the central theory directly: shuffle DHS sequences internally. Preserves per-seq
composition; destroys all motifs and positional order. If shuffled DHS ≈ unshuffled
DHS, composition is everything the model learns from biology.
