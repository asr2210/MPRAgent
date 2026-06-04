"""
015_pure_uniform_random
=======================
50k 200bp purely uniform random sequences (no genome at all). Critical
baseline: if this scores ~0.045 like random_genomic, the eval is largely
INSENSITIVE to genome content — it just rewards "diverse, well-distributed
50k×200bp libraries". If it crashes < 0.03, genome matters.
"""
import sys
from pathlib import Path

import numpy as np

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200


def main():
    rng = np.random.default_rng(SEED)
    bases = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(N_SEQS, SEQ_LEN))
    arr = bases[idx]
    seqs = ["".join(row) for row in arr]
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
