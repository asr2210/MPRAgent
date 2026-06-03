"""
005_satmut
==========
Saturation mutagenesis-style library: 100 "seed" sequences from cCREs,
each replicated as 500 variants differing from the seed by 1-3 random
single-base substitutions. Total: 50,000 sequences.

This is a fundamentally different design from random + biological libraries.
The model trained on this learns position-by-position activity effects within
a fixed set of contexts. If the eval rewards "what does a single mutation do",
this should score very differently from previous experiments.

Hypothesis: if eval is about regulatory grammar in genomic context, this
scores moderately (model learns 100 contexts well, transfers little). If
eval is about per-base position effects (saturation mutagenesis-style),
this scores much higher than other libraries.
"""
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_SEEDS = 100
VARIANTS_PER = 500
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CCRE_PATH = DATA_DIR / "ccre_v4.bed"
HG38_PATH = DATA_DIR / "hg38.fa"
CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_ccres():
    recs = []
    with open(CCRE_PATH) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            if chrom not in CANONICAL:
                continue
            start, end = int(parts[1]), int(parts[2])
            mid = (start + end) // 2
            recs.append((chrom, mid))
    return recs


def mutate(seq: str, n_muts: int, rng: np.random.Generator) -> str:
    arr = list(seq)
    L = len(arr)
    bases = ["A", "C", "G", "T"]
    positions = rng.choice(L, size=n_muts, replace=False)
    for p in positions:
        current = arr[p]
        alts = [b for b in bases if b != current]
        arr[p] = alts[int(rng.integers(0, 3))]
    return "".join(arr)


def main():
    rng = np.random.default_rng(SEED)
    print("Loading cCREs", file=sys.stderr)
    recs = load_ccres()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    # Pick N_SEEDS valid seeds (sequence with no N).
    seeds = []
    perm = rng.permutation(len(recs))
    pi = 0
    while len(seeds) < N_SEEDS and pi < len(perm):
        chrom, mid = recs[perm[pi]]
        pi += 1
        L = len(fa[chrom])
        start = max(0, mid - HALF)
        end = start + SEQ_LEN
        if end > L:
            end = L
            start = end - SEQ_LEN
        if start < 0:
            continue
        s = str(fa[chrom][start:end])
        if len(s) != SEQ_LEN or set(s) - set("ACGT"):
            continue
        seeds.append(s)
    print(f"Selected {len(seeds)} seeds", file=sys.stderr)

    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for seed in seeds:
            for _ in range(VARIANTS_PER):
                n_muts = int(rng.integers(1, 4))  # 1, 2, or 3 mutations
                v = mutate(seed, n_muts, rng)
                f.write(v + "\n")
    print(f"Wrote {N_SEEDS * VARIANTS_PER} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
