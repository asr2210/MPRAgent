"""
Experiment 027: Multi-seed Dirichlet(0.5) union (10 seeds x 5k each).

Combines 5k samples from each of seeds 0, 1, 2, ..., 9 = 50k total.
Tests whether averaging across seeds smooths out luck (giving ~mean 0.137)
or whether the BEST individual compositions from each seed dominate
(potentially > mean of individual seeds).

Theoretically equivalent to a single 50k Dirichlet(0.5) sample if all
seeds are i.i.d. — should give ~0.1371 ± 0.001 (very low variance).

If 027 substantially > 0.1371: union of seeds produces a different
composition distribution that the evaluator prefers.
"""
from pathlib import Path
import numpy as np

OUT = Path(__file__).resolve().parent / "sequences_0.txt"
N_PER_SEED = 5_000
SEEDS = list(range(10))
LEN = 200
BASES = np.array(list("ACGT"))


def main():
    all_seqs = []
    for s in SEEDS:
        rng = np.random.default_rng(s)
        probs = rng.dirichlet((0.5,)*4, size=N_PER_SEED)
        for i in range(N_PER_SEED):
            all_seqs.append("".join(BASES[rng.choice(4, size=LEN, p=probs[i])]))
    rng_shuf = np.random.default_rng(42)
    rng_shuf.shuffle(all_seqs)
    OUT.write_text("\n".join(all_seqs) + "\n")
    print(f"wrote {len(all_seqs)} multi-seed dirichlet seqs to {OUT}")


if __name__ == "__main__":
    main()
