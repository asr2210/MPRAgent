# Experiment 011 — Hyperactive design via simulated annealing

## Design
50k parallel SA trajectories, each starting from a random GC-50 200bp
scaffold. At each of 200 steps, propose a single-base mutation, accept
with probability  exp(delta / T)  where delta = Δ Malinois max(K562, HepG2,
SKNSH). Cooling schedule T: 1.0 → 0.05 linear.

Selection stats:
- Mean final Malinois score: 5.51 (vs random ~0.7, 007 top 2.60, 011 max 15.10)
- GC: mean 0.510, std 0.034 — stayed near 0.50 (SA naturally avoids extreme GC
  because Malinois drops at extreme composition).

## Hypothesis
Sequences predicted FAR beyond what random sampling can find should contain
dense, optimally-spaced TF binding sites. Training on these should expose
the model to the GRAMMAR of activation, transferable across cell types.

## Result
- eval_01 = **0.3947** (vs 007 = 0.3969, gc_50 = 0.397)
- Mean across 14 evals = **0.3839** (vs 007 = 0.3861)
- Per-cell-type on eval_01: K562 0.617, HepG2 0.432, SK-N-SH 0.136
- Runtime: 1200s

## Interpretation
**Hyperactive sequences are NOT a new lever.** Despite 2× higher predicted
activity than 007's top-magnitude selection, the trained model scored
slightly WORSE.

The pattern: there appears to be an OPTIMAL training activity level around
mean Malinois 2-3 (007). Going higher (cCRE+oracle 4.1, hyperactive 5.5)
makes the model worse — likely because the Malinois optimization landscape
diverges from real-MPRA landscape at extreme predictions.

Equivalently: SA finds sequences with stereotypical TF binding clusters
that Malinois rewards heavily but that may be UNREALISTIC ARRANGEMENTS not
seen in MPRA test sequences. The trained model overfits to these unrealistic
patterns.

## What this rules out
- Hyperactive design as a way to push past 0.397
- The "more extreme = better training signal" hypothesis
- Both ends of the activity spectrum hurt vs moderate selection

## Theory v6 reinforced
The 0.397 ceiling is robust. Single-strategy training-data selection
saturates here, regardless of whether the strategy is:
- Random
- Real cCREs
- Motif-implanted
- Oracle-selected (any criterion)
- Hyperactive designed

To break, we likely need either:
1. Mixed-strategy diversity (combining multiple distinct subpools)
2. Acceptance that ceiling is evaluator-structural and shift focus to
   stress-testing the limits of what each strategy explains.

Best library remains 007 random+Malinois top: eval_01 = 0.3969, mean_r = 0.3861.
