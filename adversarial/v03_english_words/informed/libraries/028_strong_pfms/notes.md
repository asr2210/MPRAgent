# 028 — Strong PFMs only (top 30% by info content)

## Plan
Filter 289 → 86 PFMs by mean information content (entropy filter).
Keep crispest, most-unambiguous motifs. Use exp 020 design: 8 motifs
overlap, uniform random backbone, 50k seqs.

## Result
**eval_01 = 0.4205.** K562 0.588, HepG2 0.618, SK-N-SH 0.055.
14-eval avg: 0.4205. Below 020 plateau (0.4284).

## What this teaches
- Filtering weak PFMs slightly HURTS (−0.008). Low-IC PFMs are not
  diluting signal — they're contributing useful diversity.
- Diversity > "quality" of individual PFMs.
- The 86 strongest PFMs are heavily biased toward AP-1/MYB/MAX style
  short, GC-rich motifs (top 5: OLIG2, FOS-1, MYB, MAX, MYB23). This
  family-skew may explain the slight K562/HepG2 boost and SK-N-SH
  drop — narrower repertoire = worse generalization.

## Pattern: ALL pool-restriction strategies underperform 289-full
- 011 (consensus): 0.4169
- 012 (same-tf clustered): 0.4196
- 013 (focused pool 142): 0.4220
- 014 (clean exact pool): 0.4184
- 028 (strong PFMs 86): 0.4205

Versus 008 (289 full, density 5 nooverlap): 0.4283
Versus 020 (289 full, density 8 overlap): 0.4284

**Diversity within the 289 pool is contributing real, non-redundant
signal.** The plateau isn't from pool noise.

## Next
Exp 029: candidate-pool diversity selection. Generate 150k seqs from
the BEST known designs (mix 020-style and 026-style), score by 6-mer
coverage (proxy for diversity), keep most diverse 50k. Pure data-side
approach — bypasses the per-cell-type opt-conflict entirely by letting
the library composition emerge from a diversity criterion.
