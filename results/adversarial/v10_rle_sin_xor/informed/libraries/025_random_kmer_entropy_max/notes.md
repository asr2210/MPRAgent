# Exp 025: top-50k 6-mer entropy from 500k random
- eval_01=**0.5180**  K562=0.988  HepG2=0.5617  SKNSH=0.0042

Selecting for k-mer entropy from random pool slightly REDUCES K562 (0.995→0.988).
Any non-uniform selection from a random pool moves K562 away from saturation.
HepG2 unchanged. SK-N-SH marginally positive (lucky).

Insight: even "implicitly biological" selection (high entropy) costs K562 points.
The K562 oracle truly wants uniform random selection.
