"""
Experiment 015: Block-Dirichlet.

Each 200bp sequence = 4 blocks of 50bp. Each block has its own
Dirichlet(0.5) composition. Bases iid within each block.

Tests discrete within-sequence composition variation. Differs from:
- 003 (Markov-Dirichlet, 0.135): no Markov transitions, just block changes
- 009 (gradient, 0.137): blocks are sharp not smooth

If model truly is bag-of-bases, block-Dirichlet matches pure Dirichlet
(0.138 ± 0.003 noise). If it cares about local structure, blocks
either help (added diversity) or hurt (broken local context).

Last clean test of within-sequence structure axis before moving to
seed-scanning and final refinement.
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
BLOCKS = 4
BLOCK_LEN = LEN // BLOCKS  # 50
ALPHA = (0.5, 0.5, 0.5, 0.5)
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    probs = rng.dirichlet(ALPHA, size=(N, BLOCKS))  # (N, B, 4)
    seqs = []
    for i in range(N):
        parts = []
        for j in range(BLOCKS):
            b = rng.choice(4, size=BLOCK_LEN, p=probs[i, j])
            parts.append("".join(BASES[b]))
        seqs.append("".join(parts))
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} block-dirichlet seqs (4x50bp blocks) to {OUT}")


if __name__ == "__main__":
    main()
