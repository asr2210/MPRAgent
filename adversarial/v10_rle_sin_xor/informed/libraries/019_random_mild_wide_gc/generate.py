"""
Experiment 019: mild widening of GC: uniform GC ∈ [0.4, 0.6] per sequence.

Standard random has GC ~0.5 ± 0.035. This widens slightly to ~0.5 ± 0.06 uniform.
Tests if a small increase in compositional variance is positive.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

rng = np.random.default_rng(SEED)
with open(OUT, "w") as f:
    for _ in range(50_000):
        gc_target = rng.uniform(0.4, 0.6)
        is_gc = rng.random(200) < gc_target
        seq = []
        for ig in is_gc:
            if ig:
                seq.append("G" if rng.random() < 0.5 else "C")
            else:
                seq.append("A" if rng.random() < 0.5 else "T")
        f.write("".join(seq) + "\n")
print("Done")
