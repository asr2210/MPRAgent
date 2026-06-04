"""
E5: gc_50 strict replication.

For each of 50,000 seqs: randomly pick 100 positions out of 200 to be G/C,
remaining 100 to be A/T. Within each GC position, G or C with P=0.5 (independent).
Within each AT position, A or T with P=0.5 (independent).

Replicates published gc_50 baseline (5-seed mean: 0.8591 eval_01).
Tests whether tightening per-seq GC to exactly 100 (while keeping G:C and
A:T free) lifts marginally over random_uniform.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
N_GC = 100  # exactly 100 G/C bases

rng = np.random.default_rng(SEED)

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    with open(out_path, "w") as f:
        for _ in range(N_SEQS):
            seq = np.empty(SEQ_LEN, dtype="<U1")
            # pick 100 random positions to be GC
            gc_pos = rng.choice(SEQ_LEN, size=N_GC, replace=False)
            at_mask = np.ones(SEQ_LEN, dtype=bool)
            at_mask[gc_pos] = False
            # assign G or C to gc_pos
            gc_bases = rng.choice(["G", "C"], size=N_GC)
            seq[gc_pos] = gc_bases
            # assign A or T to at_pos
            at_bases = rng.choice(["A", "T"], size=SEQ_LEN - N_GC)
            seq[at_mask] = at_bases
            f.write("".join(seq) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
