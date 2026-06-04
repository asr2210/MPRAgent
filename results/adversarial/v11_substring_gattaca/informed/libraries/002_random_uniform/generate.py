"""
E2: random_uniform — baseline calibration.

50,000 IID random 200bp sequences. P(A)=P(C)=P(G)=P(T)=0.25.
Single random seed (0). Confirms my pipeline matches the published
random_uniform baseline (~0.857 mean_r at eval_01).
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200

rng = np.random.default_rng(SEED)

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    # Sample all at once: shape (N_SEQS, SEQ_LEN)
    arr = rng.integers(0, 4, size=(N_SEQS, SEQ_LEN), dtype=np.uint8)
    bases = np.array(["A", "C", "G", "T"])
    with open(out_path, "w") as f:
        for row in arr:
            f.write("".join(bases[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
