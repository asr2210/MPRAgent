# Experiment 004 — asymmetric cCRE + PLS boost + paired flanks

## Design
- Positives: 10k dELS, 4k pELS, 3k PLS (BOOST), 4k CTCF-only, 4k DNH3
- 25k paired flanks (±1500-3000bp from centers)
- Hypothesis: PLS boost helps eval_13 + cross-cell-type generalization
  (per LentiMPRA: promoters share TFs across cell types)

## Result — eval_01 = 0.1479
- K562_r = -0.003 / SKNSH = 0.450
- Slightly better than exp 003 (balanced) = 0.145
- But far from blind 013 (asymm + no PLS) = 0.173

## What this tells me
PLS boost did not pay off as expected.
- eval_13 = 0.144 in 004 vs 0.126 in blind 013 — a small win, but well
  below blind 008/011's eval_13 ≈ 0.158
- The extra PLS slots (3k vs blind's ~0.6k) cost ~2.4k dELS/CTCF/DNH3
  slots, which seems to have hurt overall

Effective composition contrast (blind 013 vs my 004):
  | class | blind 013 | my 004 |
  |-------|-----------|--------|
  | dELS  | ~11k      | 10k    |
  | pELS  | ~2.4k     | 4k     |
  | PLS   | ~0.6k     | 3k     |
  | CTCF  | ~5.5k     | 4k     |
  | DNH3  | ~5.35k    | 4k     |

The bulk of blind's win comes from CTCF/DNH3 boost. My downweighting of
those was suboptimal.

## Theory v4 update
PLS isn't the magic ingredient for eval_13. eval_13 may be measuring
something else (a specific tissue / cell line not in K562/HepG2/SKNSH?
A particular motif class?). Need to investigate eval-set patterns more
carefully.

The CTCF/DNH3 boost in blind 013 is the load-bearing ingredient. Those are
rare in nature but functionally distinct (CTCF binding = insulator; DNH3 =
proximal regulatory). Each contributes complementary motif signal.

## Next direction
I've spent 4 experiments on cCRE+flank variations. Time to try something
fundamentally different. Options:
- Use the published Malinois MPRA dataset (Table_S2, ~798k sequences with
  measured K562/HepG2/SKNSH activity) — direct empirical signal
- Add synthetic motif-loaded sequences (JASPAR-driven)
- Active learning: train, score, select uncertain sequences
