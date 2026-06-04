"""
Exp 029: mixture of 10k each from seeds {1, 7, 42, 100, 2024}.

If "random" libraries have any seed-specific bias, mixing should produce a
median-like result. If results are seed-independent, mixture ≈ any individual.
Indirect test for whether single-seed variance is genuine RNG luck vs systematic.
"""
import os, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")

bases = np.array(list("ACGT"))
seqs = []
for s in [1, 7, 42, 100, 2024]:
    rng = np.random.default_rng(s)
    for _ in range(10_000):
        seqs.append("".join(bases[rng.integers(0, 4, 200)]))

rng = np.random.default_rng(0)
rng.shuffle(seqs)
assert len(seqs) == 50_000
with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print("Done")
