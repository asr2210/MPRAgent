"""
Experiment 013: Trimmed Dirichlet(0.5).

Sample Dirichlet(0.5), reject any composition where min(p) < 0.03.
This removes the near-homopolymer tail (sequences where one base has
<3% representation, i.e. <6 occurrences in 200bp). Keeps most of the
distribution intact.

Hypothesis: the very-skewed tail of Dirichlet(0.5) is the unproductive
part. Trimming should:
- Keep the wide composition spread that helps the model
- Avoid the near-homopolymer sequences (we know those hurt:
  v06 homopolymer_rich = 0.058)
- Slightly tighten the distribution toward "moderately diverse"

Reject ratio with min_p=0.03 is ~15-20% empirically with alpha=0.5,
so we need ~65k candidates to get 50k accepted.
"""
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
ALPHA = (0.5, 0.5, 0.5, 0.5)
MIN_P = 0.03
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)

    # Oversample then trim
    over = 200_000
    cand = rng.dirichlet(ALPHA, size=over)
    keep = cand.min(axis=1) >= MIN_P
    print(f"acceptance rate: {keep.mean():.3f}")
    cand = cand[keep]
    assert len(cand) >= N, f"not enough accepted: {len(cand)}"
    probs = cand[:N]

    seqs = []
    for i in range(N):
        b = rng.choice(4, size=LEN, p=probs[i])
        seqs.append("".join(BASES[b]))
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} trimmed-dirichlet seqs to {OUT}")
    print(f"composition stats: GC mean={probs[:,[1,2]].sum(1).mean():.3f} "
          f"std={probs[:,[1,2]].sum(1).std():.3f}")
    print(f"max single base mean={probs.max(axis=1).mean():.3f}")


if __name__ == "__main__":
    main()
