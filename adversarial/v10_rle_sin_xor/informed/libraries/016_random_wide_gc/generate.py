"""
Experiment 016: random with WIDE per-sequence GC variance (uniform [0.2, 0.8]).

Following major exp 015 finding: HepG2 is composition-variance driven.
This experiment AMPLIFIES per-sequence GC variance vs natural random's ~0.035 stdev.
Each sequence gets a target GC drawn uniformly from [0.2, 0.8], then bases sampled
to match that target.

Hypothesis: more compositional variance → more HepG2 signal → higher mean_r.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

rng = np.random.default_rng(SEED)
with open(OUT, "w") as f:
    for _ in range(50_000):
        gc_target = rng.uniform(0.2, 0.8)
        # Each position: P(G or C) = gc_target; then equal split GC vs AT
        is_gc = rng.random(200) < gc_target
        # Within GC pick G/C, within AT pick A/T
        seq = []
        for ig in is_gc:
            if ig:
                seq.append("G" if rng.random() < 0.5 else "C")
            else:
                seq.append("A" if rng.random() < 0.5 else "T")
        f.write("".join(seq) + "\n")
print("Done")
