"""
Experiment 011: random sequences from human-genome dinucleotide Markov-1 chain.

Hypothesis: maybe pure i.i.d. random misses natural local context. Markov-1
sequences have natural dinucleotide frequencies (CpG suppression, etc.) but
otherwise random. Tests if natural-feel local statistics help.

If better than i.i.d. random → natural dinucleotide context matters.
If worse → uniform is optimal.
"""
import os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1

# Approximate human-genome dinucleotide transitions (rough from literature)
# rows = current base, cols = next base. CpG dramatically suppressed.
# Source: rough human-genome dinucleotide frequencies.
trans = np.array([
    # to: A      C      G      T
    [0.298, 0.198, 0.202, 0.302],  # from A
    [0.328, 0.218, 0.022, 0.432],  # from C (CpG suppressed)
    [0.292, 0.241, 0.222, 0.245],  # from G
    [0.215, 0.211, 0.288, 0.286],  # from T
])
# Normalize
trans = trans / trans.sum(axis=1, keepdims=True)

# Initial distribution: rough human genome
init = np.array([0.295, 0.205, 0.205, 0.295])
init /= init.sum()

bases_idx = "ACGT"
rng = np.random.default_rng(SEED)
with open(OUT, "w") as f:
    for _ in range(50_000):
        seq = []
        cur = rng.choice(4, p=init)
        seq.append(bases_idx[cur])
        for _ in range(199):
            cur = rng.choice(4, p=trans[cur])
            seq.append(bases_idx[cur])
        f.write("".join(seq) + "\n")
print("Done")
