"""
Experiment 025: Dirichlet(0.5) + small infusion of K562-relevant synthseqs.

45k Dirichlet(0.5) + 5k synthseqs from component 15 (Myeloid/erythroid,
K562-relevant). Tests if a small dose of real K562-relevant sequences
boosts K562 head (~0.04 → higher) without sinking HepG2/SK-N-SH.

Component 15 has ~10k synthseqs; sample 5k of them.
"""
from pathlib import Path
import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "data"
OUT = Path(__file__).resolve().parent / "sequences_0.txt"
SEED = 42
N = 50_000
N_DIR = 45_000
N_REAL = 5_000
LEN = 200
BASES = np.array(list("ACGT"))


def main():
    rng = np.random.default_rng(SEED)

    syn = pd.read_csv(DATA / "train_synthseqs.csv.gz", sep="\t",
                      compression="gzip")
    pool = syn[syn["component"] == 15]
    real_idx = rng.choice(len(pool), size=N_REAL, replace=False)
    real_seqs = pool["raw_sequence"].iloc[real_idx].tolist()

    probs = rng.dirichlet((0.5,)*4, size=N_DIR)
    dir_seqs = ["".join(BASES[rng.choice(4, size=LEN, p=probs[i])])
                for i in range(N_DIR)]

    seqs = dir_seqs + real_seqs
    rng.shuffle(seqs)
    OUT.write_text("\n".join(seqs) + "\n")
    print(f"wrote {N} ({N_DIR} Dirichlet + {N_REAL} K562 synthseqs)")


if __name__ == "__main__":
    main()
