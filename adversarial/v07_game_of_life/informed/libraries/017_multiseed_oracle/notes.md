# Experiment 017 — Multi-seed oracle pool (5 × 10k)

## Design
5 independent oracle runs (SEEDS 1-5), each:
- 100k GC-50 random pool
- Score with Malinois
- Top 10k by max(K562, HepG2, SKNSH)

Concatenate the 50k into a single library.

## Result
- eval_01 = **0.3977** (Δ vs 007 = +0.0008, vs 014 = +0.0041)
- mean_r = **0.3862** (Δ vs 007 = +0.0001, vs 014 = +0.0035)
- Per-cell on eval_01: K562=0.618, HepG2=0.436, SK-N-SH=0.139

## Interpretation
**Multi-seed oracle reliably reaches 007's level — confirms 007 wasn't luck.**

The single-seed oracle distribution spans:
- 014 (unlucky seed=2): 0.3827
- 007 (lucky seed=1): 0.3861
- 017 (averaging 5 seeds): 0.3862

So the TRUE oracle-selection ceiling is mean_r ≈ 0.386 (eval_01 ≈ 0.397).
The 014 unlucky case sits ~0.003 below this. Multi-seed converges to the
upper end of the noise distribution — it provides RELIABILITY but does
NOT break the ceiling.

This is a strong signal that the ceiling is the (model+50k+evaluator)
limit when applied to GC-50 + oracle-selected sequences. No amount of
seed-averaging will exceed it.

## Theory v11
The oracle-selection ceiling is at:
- mean_r ~0.386
- eval_01 ~0.397

This ceiling is robust:
- Pool size (500k vs 4M): doesn't help (013 slightly worse)
- Selection criterion (max vs span vs discrim): all converge at noise level
- Seed averaging: converges to ceiling reliably but doesn't exceed it

To break it, we need a fundamentally different sequence-level information
content. Options:
1. **Pan-active oracle** (018): select by MIN cross-cell — bias for
   sequences that activate ALL 3 cells (presumably more transferable).
2. **Cross-model oracle**: use Malinois + a different sequence model
   for agreement-based selection.
3. **Real natural distribution**: cCREs with no GC filter (natural
   composition variation per sequence within natural-like subspace).
4. **Sequence subspace expansion**: add reverse-complements to expand
   effective training set within the same library size.

## What this rules out
- Single-seed luck explaining 007's higher result
- Multi-seed averaging as a way to break the ceiling
- Hope that more replicates within the same recipe would help
