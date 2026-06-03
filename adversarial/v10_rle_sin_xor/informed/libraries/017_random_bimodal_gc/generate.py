"""
Experiment 017: bimodal GC distribution: half sequences GC=0.3, half GC=0.7.

If HepG2 cares about composition variance specifically (not the variance shape),
bimodal might provide cleaner signal: two clear classes the oracle can separate.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

rng = np.random.default_rng(SEED)
with open(OUT, "w") as f:
    for i in range(50_000):
        gc_target = 0.3 if i < 25_000 else 0.7
        is_gc = rng.random(200) < gc_target
        seq = []
        for ig in is_gc:
            if ig:
                seq.append("G" if rng.random() < 0.5 else "C")
            else:
                seq.append("A" if rng.random() < 0.5 else "T")
        f.write("".join(seq) + "\n")
print("Done")
