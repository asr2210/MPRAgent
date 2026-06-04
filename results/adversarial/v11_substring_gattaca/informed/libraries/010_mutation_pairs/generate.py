"""
E10: mutation_pairs

25,000 IID random uniform sequences (parents) + 25,000 single-point
mutations (children). Each child = parent with one position changed to
a random different base. Tests contrastive pair augmentation.

Per-position distribution remains IID uniform.
"""
import numpy as np
import os

SEED = 0
N_PAIRS = 25_000
N_SEQS = 50_000
SEQ_LEN = 200

rng = np.random.default_rng(SEED)

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    parents = rng.integers(0, 4, size=(N_PAIRS, SEQ_LEN), dtype=np.int8)
    children = parents.copy()
    mut_pos = rng.integers(0, SEQ_LEN, size=N_PAIRS)
    for n in range(N_PAIRS):
        old = parents[n, mut_pos[n]]
        # pick a new base not equal to old
        choices = [b for b in range(4) if b != old]
        new = rng.choice(choices)
        children[n, mut_pos[n]] = new
    full = np.vstack([parents, children])
    perm = rng.permutation(N_SEQS)
    full = full[perm]
    bases = np.array(["A", "C", "G", "T"])
    with open(out_path, "w") as f:
        for row in full:
            f.write("".join(bases[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
