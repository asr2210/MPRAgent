# 004 — Dirichlet(0.3) extreme composition

## What I built
50,000 sequences, each with base frequencies ~ Dirichlet(0.3,0.3,0.3,0.3) (concentrates near simplex corners → many sequences dominated by 1-2 bases), 200bp i.i.d. from that composition.

## Predicted
If composition variance is the lever (H8 from theory v1), more extreme per-sequence compositions should beat Dirichlet(1.0) (baseline 0.0768) on eval_01.

## Result
- **eval_01 = 0.0786** — beats Dirichlet(1.0) baseline (0.0768) and all my previous experiments. NEW BEST.
- eval_08 = 0.0716 (vs 0.0728 in exp 002, 0.0627 in exp 001) — close to best
- mean across evals ≈ 0.0954 — best so far
- Cell-type breakdown roughly K562 ≈ HepG2 > SKNSH consistently

## Theory updates
- T1.2 (composition variance lever) confirmed by 0.0786 > 0.0768.
- T1.4 (~0.08 ceiling) close but exp 004 nudges past 0.0786 — small room remains.

## Next experiment
EXP 005 — Dirichlet(0.1) extreme. Tests monotonicity of the alpha → score relationship. Two possible outcomes:
- (a) 0.1 > 0.3 > 1.0: monotonic. Push toward more extreme compositions.
- (b) 0.1 < 0.3 < 1.0 — implies sweet spot. Optimal alpha somewhere in (0.1, 1.0).

If (b), homopolymer_rich (0.0570) hints that pure-homopolymers hurt. So there must be a balance.

Predicted: probably (b) — Dirichlet(0.1) puts so much mass at corners that many sequences are near-homopolymer, which strategies.md showed hurt.
