"""
Experiment 013: 25k random uniform + their 25k reverse complements.

Hypothesis: The model may benefit from strand symmetry — seeing both strands
of each random sequence teaches strand-invariant pattern matching with the
same total entropy as fully random 50k. If revcomp augmentation helps,
mean_r should equal or exceed pure random.

If equal: 50k unique vs 25k+25k revcomp is information-equivalent.
If worse: strand redundancy hurts diversity.
If better: strand symmetry is a useful inductive bias.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

rng = np.random.default_rng(SEED)
bases = np.array(list("ACGT"))
COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}

seqs = []
for _ in range(25_000):
    s = "".join(bases[rng.integers(0, 4, 200)])
    seqs.append(s)
    seqs.append("".join(COMP[b] for b in s[::-1]))

assert len(seqs) == 50_000
rng.shuffle(seqs)
with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"Wrote {len(seqs)} sequences (25k random + 25k revcomp)")
