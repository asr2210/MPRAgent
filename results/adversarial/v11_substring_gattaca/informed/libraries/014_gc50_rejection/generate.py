"""
E14: gc50_rejection

Pure IID random rejection-sampled to per-seq GC count in {99, 100, 101}.
Unlike E7 (which flipped bases to force exact GC=100, breaking per-position
marginals), this keeps only sequences whose GC happens to land tight.

Tests if the small published gc_50 vs random_uniform edge (+0.0025) is
reproducible signal or noise.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 99, 101  # very tight band around 100 (50%)

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    # Acceptance rate: P(GC in {99,100,101}) ~ 17% for Binomial(200, 0.5)
    # Pool 500k to be safe.
    BATCH = 500_000
    selected = []
    total_drawn = 0
    while len(selected) < N_SEQS:
        batch = rng.integers(0, 4, size=(BATCH, SEQ_LEN), dtype=np.int8)
        total_drawn += BATCH
        # G=2, C=1 in BASES indexing... but indices: A=0,C=1,G=2,T=3
        # GC count = sum(batch == 1) + sum(batch == 2) per row
        is_gc = (batch == 1) | (batch == 2)
        gc_counts = is_gc.sum(axis=1)
        mask = (gc_counts >= GC_LOW) & (gc_counts <= GC_HIGH)
        kept = batch[mask]
        selected.append(kept)
        total_kept = sum(len(s) for s in selected)
        print(f"  Drew {total_drawn}, kept {total_kept} (target {N_SEQS}). Acceptance ~ {total_kept/total_drawn:.3f}")
        if total_kept >= N_SEQS:
            break
    seqs = np.vstack(selected)[:N_SEQS]
    print(f"Final: {len(seqs)} sequences from {total_drawn} drawn ({len(seqs)/total_drawn:.3f} acceptance).")
    # Verify
    is_gc = (seqs == 1) | (seqs == 2)
    gc = is_gc.sum(axis=1)
    print(f"GC count: min={gc.min()}, max={gc.max()}, mean={gc.mean():.3f}")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
