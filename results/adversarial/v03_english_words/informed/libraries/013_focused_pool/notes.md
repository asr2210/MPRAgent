# 013 — Focused 50-motif pool (canonical cell-type TFs)

## Plan
Drop universals and rare TFs; keep only the 31 canonical TF tokens for
K562/HepG2/SK-N-SH. 3 motifs/seq in random uniform backbone.

## Result
**eval_01 = 0.4220.** Worse than exp 008 (0.4283) by -0.006.
K562: 0.591, HepG2: 0.621, SK-N-SH: 0.055 (slightly above floor).

## What went wrong
The substring token list pulled in a LOT of plant MYBs (MYB1, MYB10...MYB99,
MYB3R1, etc. — Arabidopsis transcription factors) and extended KLF/SOX
family members not relevant to K562/HepG2/SK-N-SH. Final pool was 142 PFMs,
not 50 — but a meaningfully different mix than exp 008. Many are plant TFs
binding random-looking motifs.

So this experiment doesn't cleanly test "smaller focused pool" — it tests
"different (partly polluted) pool of similar size." It still lost to 008,
which is informative: arbitrary motif identity matters less than coverage
of human/mammalian regulators.

SK-N-SH bumped to 0.055 (highest of all experiments). Possibly noise, but
worth noting that of the 14 evaluations, SK-N-SH peaked here. This is the
first hint SK-N-SH might be moveable.

## Theory T10 (refined)
- Pool composition matters. 289 mixed human regulators (exp 008) beats 142
  with plant-TF pollution.
- Universals are NOT noise — removing them and adding plant MYBs lost
  ground. The 008 pool's mix of cell-type + universals was useful.
- T10 still holds: TF presence is the surrogate's feature, and more
  human-relevant TFs per sequence is better.

## Next
- Exp 014: CLEAN small pool. Use exact name matching to get only the
  intended ~30 human TFs (no plant MYBs). Tests whether smaller-but-clean
  beats larger-but-mixed.
- Or pivot: exp 014 = mixture (motif + random uniform) to add sequence
  variety.
