"""
E26: per_position_stratified

For each position, create a column with EXACTLY 12500 A, 12500 C,
12500 G, 12500 T shuffled. Stack 200 columns to form 50k sequences.

No GC filter applied. Per-position is exactly uniform; per-seq GC is
Binomial(200, 0.5) — same as random_uniform.

Tests whether per-position EXACT balance alone is the win, vs the GC
filter in E15.
"""
import numpy as np
import os

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
TARGET = N_SEQS // 4  # 12500 of each base

rng = np.random.default_rng(SEED)
BASES = np.array(["A", "C", "G", "T"])

def main():
    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    print(f"Building {SEQ_LEN} stratified columns of {N_SEQS} bases each ({TARGET} each base)...")
    seqs = np.empty((N_SEQS, SEQ_LEN), dtype=np.int8)
    base_template = np.concatenate([
        np.full(TARGET, 0, dtype=np.int8),
        np.full(TARGET, 1, dtype=np.int8),
        np.full(TARGET, 2, dtype=np.int8),
        np.full(TARGET, 3, dtype=np.int8),
    ])
    for p in range(SEQ_LEN):
        col = base_template.copy()
        rng.shuffle(col)
        seqs[:, p] = col
    # Verify
    is_gc = (seqs == 1) | (seqs == 2)
    gc = is_gc.sum(axis=1)
    print(f"GC: min={gc.min()}, max={gc.max()}, mean={gc.mean():.2f}, sd={gc.std():.2f}")
    # Per-position counts
    for p in [0, 100, 199]:
        col = seqs[:, p]
        counts = [int((col == b).sum()) for b in range(4)]
        print(f"  pos {p}: A={counts[0]}, C={counts[1]}, G={counts[2]}, T={counts[3]}")
    rng.shuffle(seqs)
    with open(out_path, "w") as f:
        for row in seqs:
            f.write("".join(BASES[row]) + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
