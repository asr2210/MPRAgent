"""
E7: true_gc50 — IID per-position sampling with rebalance to exact GC=100.

For each of 50k seqs:
1. Sample each position IID from {A,C,G,T} with P=0.25.
2. If GC count != 100, randomly flip excess GC positions to A/T (or vice versa).

Preserves per-position randomness (which SK-N-SH evidently needs) while
enforcing exact per-seq GC=100 (the marginal lift from published gc_50).
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200

rng = np.random.default_rng(SEED)

BASES = np.array(["A", "C", "G", "T"])
GC = {"G", "C"}
AT = {"A", "T"}

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    arr = rng.integers(0, 4, size=(N_SEQS, SEQ_LEN), dtype=np.int8)
    # 0=A, 1=C, 2=G, 3=T. GC mask = (arr==1) | (arr==2)
    is_gc = (arr == 1) | (arr == 2)
    gc_counts = is_gc.sum(axis=1)
    with open(out_path, "w") as f:
        for n in range(N_SEQS):
            seq = arr[n].copy()
            gc_count = int(gc_counts[n])
            if gc_count > 100:
                # flip (gc_count - 100) random GC positions to AT
                gc_idx = np.where((seq == 1) | (seq == 2))[0]
                to_flip = rng.choice(gc_idx, size=gc_count - 100, replace=False)
                new_bases = rng.choice([0, 3], size=len(to_flip))  # A or T
                seq[to_flip] = new_bases
            elif gc_count < 100:
                at_idx = np.where((seq == 0) | (seq == 3))[0]
                to_flip = rng.choice(at_idx, size=100 - gc_count, replace=False)
                new_bases = rng.choice([1, 2], size=len(to_flip))  # C or G
                seq[to_flip] = new_bases
            # Verify
            assert ((seq == 1) | (seq == 2)).sum() == 100
            f.write("".join(BASES[seq]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
