"""
E13: low_kmer_entropy

Symmetric opposite of E12. Generate 200k random uniform sequences,
select BOTTOM 50k by within-sequence 6-mer diversity (most repetitive).
If E13 hurts where E12 was flat, the eval detects micro-distribution
mismatch. If E13 also stays at ~0.84, the eval is purely macro-distributional.
"""
import numpy as np
import os

SEED = 0
N_POOL = 200_000
N_SEQS = 50_000
SEQ_LEN = 200
K = 6

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    print(f"Generating {N_POOL} random sequences...")
    pool = rng.integers(0, 4, size=(N_POOL, SEQ_LEN), dtype=np.int8)
    print("Computing 6-mer diversity per sequence...")
    powers = 4 ** np.arange(K - 1, -1, -1, dtype=np.int64)
    n_kmers = SEQ_LEN - K + 1
    kmer_vals = np.zeros((N_POOL, n_kmers), dtype=np.int64)
    for i in range(K):
        kmer_vals += pool[:, i : i + n_kmers].astype(np.int64) * powers[i]
    unique_counts = np.empty(N_POOL, dtype=np.int32)
    for r in range(N_POOL):
        unique_counts[r] = len(np.unique(kmer_vals[r]))
    print(f"Unique 6-mer counts: min={unique_counts.min()}, max={unique_counts.max()}, mean={unique_counts.mean():.1f}")
    # Select BOTTOM N_SEQS (most repetitive)
    bot_idx = np.argpartition(unique_counts, N_SEQS)[:N_SEQS]
    selected = pool[bot_idx]
    print(f"Selected bottom {N_SEQS}, threshold unique 6-mers: {unique_counts[bot_idx].max()}")
    rng.shuffle(selected)
    with open(out_path, "w") as f:
        for row in selected:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
