"""
E3: balanced_bases — strict per-nucleotide balance at 50% GC.

Each sequence is a random permutation of exactly 50A + 50C + 50G + 50T.
Tests whether removing within-sequence base count variance improves over
random_uniform / gc_50.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200

rng = np.random.default_rng(SEED)

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    base_template = np.array(["A"] * 50 + ["C"] * 50 + ["G"] * 50 + ["T"] * 50)
    with open(out_path, "w") as f:
        for _ in range(N_SEQS):
            perm = rng.permutation(SEQ_LEN)
            seq = base_template[perm]
            f.write("".join(seq) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
