"""
E12: high_kmer_entropy

Generate 200k random uniform sequences. Select top 50k by within-seq
6-mer DIVERSITY (number of unique 6-mers per sequence). Tests within-IID
selection effect.
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

def count_unique_kmers(seq_int, k=6):
    """Count unique k-mers in a sequence given as int array (0-3)."""
    # Encode k-mers as base-4 integers
    powers = 4 ** np.arange(k - 1, -1, -1, dtype=np.int64)
    # Sliding window k-mer integer values
    n = len(seq_int)
    vals = np.zeros(n - k + 1, dtype=np.int64)
    for i in range(k):
        vals += seq_int[i : n - k + 1 + i].astype(np.int64) * powers[i]
    return len(np.unique(vals))

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    print(f"Generating {N_POOL} random sequences...")
    pool = rng.integers(0, 4, size=(N_POOL, SEQ_LEN), dtype=np.int8)
    print("Computing 6-mer diversity per sequence...")
    # Vectorized kmer encoding
    powers = 4 ** np.arange(K - 1, -1, -1, dtype=np.int64)
    n_kmers = SEQ_LEN - K + 1
    kmer_vals = np.zeros((N_POOL, n_kmers), dtype=np.int64)
    for i in range(K):
        kmer_vals += pool[:, i : i + n_kmers].astype(np.int64) * powers[i]
    # Compute unique count per row using sort
    unique_counts = np.empty(N_POOL, dtype=np.int32)
    for r in range(N_POOL):
        unique_counts[r] = len(np.unique(kmer_vals[r]))
    print(f"Unique 6-mer counts: min={unique_counts.min()}, max={unique_counts.max()}, mean={unique_counts.mean():.1f}")
    # Select top N_SEQS
    top_idx = np.argpartition(-unique_counts, N_SEQS)[:N_SEQS]
    selected = pool[top_idx]
    print(f"Selected {N_SEQS}, threshold unique 6-mers: {unique_counts[top_idx].min()}")
    # Shuffle selected to randomize order
    rng.shuffle(selected)
    with open(out_path, "w") as f:
        for row in selected:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
