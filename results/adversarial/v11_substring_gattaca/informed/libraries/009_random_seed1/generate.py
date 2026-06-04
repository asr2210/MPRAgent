"""
E9: random_uniform with SEED=1 — single-seed variance check.

Same as E2 but seed=1. Confirms run-to-run variance for random_uniform.
"""
import numpy as np
import os

SEED = 1
N_SEQS = 50_000
SEQ_LEN = 200

rng = np.random.default_rng(SEED)

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    arr = rng.integers(0, 4, size=(N_SEQS, SEQ_LEN), dtype=np.uint8)
    bases = np.array(["A", "C", "G", "T"])
    with open(out_path, "w") as f:
        for row in arr:
            f.write("".join(bases[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
