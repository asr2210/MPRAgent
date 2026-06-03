# E020 — Gosai chr 7, 9, 13, 19, 21, X (E008 + Malinois val chr19)

50K random Gosai from those 6 chroms (157K pool).

## Result
eval_01 = 0.3335. Per-cell K562=0.165, HepG2=0.216, SKNSH=0.620.

## Interpretation
Essentially the same as E008 (0.3359) — chr19 addition gives no boost.
Slight downward shift is just dilution noise.

So the Malinois-held-out hypothesis is partially supported (chr 7/9/13/21/X
are real test set in published Malinois; chr19 might be val and *also* helps,
but the effect is too small to detect at one seed).

The chr-lever is real but capped at +0.013, regardless of which held-out
chrom you add. Whatever signal exists is roughly saturated by E008.

## eval_04 boost
eval_04 jumped from 0.276 (E008) to 0.288 (E020). Minor diversification
benefit, but eval_01 didn't move.
