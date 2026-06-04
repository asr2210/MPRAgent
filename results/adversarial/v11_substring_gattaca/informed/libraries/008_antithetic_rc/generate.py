"""
E8: antithetic_rc

25,000 IID random uniform sequences + their 25,000 reverse complements.
Total = 50k. Preserves per-position 50% GC IID distribution exactly.

Tests whether explicit RC pair augmentation helps (if model lacks internal
RC aug) or hurts (if model has RC aug and pairs waste unique data).
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
N_BASE = N_SEQS // 2  # 25,000 base sequences, paired with their RCs
SEQ_LEN = 200

rng = np.random.default_rng(SEED)

BASES = np.array(["A", "C", "G", "T"])
COMP_INT = np.array([3, 2, 1, 0], dtype=np.int8)  # A<->T (0<->3), C<->G (1<->2)

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    arr = rng.integers(0, 4, size=(N_BASE, SEQ_LEN), dtype=np.int8)
    # Reverse complement: reverse along axis 1 then map via COMP_INT
    rc = COMP_INT[arr[:, ::-1]]
    # Stack and shuffle
    full = np.vstack([arr, rc])
    perm = rng.permutation(N_SEQS)
    full = full[perm]
    with open(out_path, "w") as f:
        for row in full:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
