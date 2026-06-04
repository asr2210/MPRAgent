"""
Exp 025: random sequences selected for maximum per-sequence 6-mer entropy.

Generate 500k random sequences. For each, compute 6-mer count entropy.
Keep top 50k by entropy. Sequences with the most uniform 6-mer distribution
have the most "information per sequence" — tests if the K562/HepG2 oracle
prefers high-entropy random samples.
"""
import os
import numpy as np
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1
N_POOL = 500_000
K = 6

rng = np.random.default_rng(SEED)
bases = np.array(list("ACGT"))

# Generate sequences and entropy scores in one pass
records = []  # (entropy, seq)
batch = 5000
n_done = 0
while n_done < N_POOL:
    bs = min(batch, N_POOL - n_done)
    arr = bases[rng.integers(0, 4, size=(bs, 200))]
    for row in arr:
        s = "".join(row)
        # 6-mer entropy: distribution over 6-mers in this sequence
        counts = Counter(s[i:i+K] for i in range(200 - K + 1))
        total = sum(counts.values())
        ps = np.array(list(counts.values()), dtype=float) / total
        H = float(-(ps * np.log2(ps)).sum())
        records.append((H, s))
    n_done += bs
print(f"Generated {len(records)} candidates")

records.sort(reverse=True)  # high entropy first
chosen = [s for _, s in records[:50_000]]
print(f"Top entropy: {records[0][0]:.4f}, 50k-th: {records[49999][0]:.4f}")

rng.shuffle(chosen)
with open(OUT, "w") as f:
    for s in chosen:
        f.write(s + "\n")
print("Done")
