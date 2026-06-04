"""
Experiment 026: Asymmetric Dirichlet (different alpha per base).

alpha = (0.4, 0.6, 0.6, 0.4) — slightly favors GC over AT, more peaky
on A and T, flatter on G and C. Mean GC fraction = 0.6 / 1.0 = 60%.

Tests whether asymmetric prior helps. Most strong baselines have GC
near 50% (gc_50 = 0.118, gc_rich at 80% = 0.058 BAD, at_rich at 20% =
0.106). So mean GC > 50% should hurt baselines. But Dirichlet's wide
spread might absorb the shift.

If 026 > 0.137: asymmetric helps — interesting!
If 026 < 0.137: symmetric Dirichlet(0.5) was optimal.
"""
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
LEN = 200
ALPHA = (0.4, 0.6, 0.6, 0.4)  # A, C, G, T
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)
    probs = rng.dirichlet(ALPHA, size=N)
    seqs = ["".join(BASES[rng.choice(4, size=LEN, p=probs[i])]) for i in range(N)]
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} dirichlet(alpha={ALPHA}) to {OUT}")
    print(f"mean composition: {probs.mean(axis=0)}")


if __name__ == "__main__":
    main()
