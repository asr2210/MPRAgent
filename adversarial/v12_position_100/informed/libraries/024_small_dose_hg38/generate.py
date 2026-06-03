#!/usr/bin/env python3
"""
Experiment 024 — 40k Dirichlet(0.3) + 10k random hg38.

Tests if small dose of biology supplements without too much dilution.
"""
import os
import gzip
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HG38 = os.path.abspath(os.path.join(HERE, "..", "..", "data", "hg38.fa.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
N_DIRICHLET = 40_000
N_BIO = 10_000
SEQ_LEN = 200
ALPHA = 0.3
ALPHABET = np.array(list("ACGT"))
USE_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX"}


def load_chroms():
    chroms = {}
    cur_name = None
    cur_chunks = []
    with gzip.open(HG38, "rt") as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if cur_name in USE_CHROMS:
                    chroms[cur_name] = "".join(cur_chunks).upper()
                cur_name = line[1:].split()[0]
                cur_chunks = []
            else:
                if cur_name in USE_CHROMS:
                    cur_chunks.append(line)
        if cur_name in USE_CHROMS and cur_name not in chroms:
            chroms[cur_name] = "".join(cur_chunks).upper()
    return chroms


def main():
    rng = np.random.default_rng(SEED)

    dirichlet_seqs = []
    for _ in range(N_DIRICHLET):
        p = rng.dirichlet([ALPHA] * 4)
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        dirichlet_seqs.append("".join(ALPHABET[idx]))
    print(f"generated {len(dirichlet_seqs)} dirichlet sequences")

    chroms = load_chroms()
    chrom_names = sorted(chroms.keys())
    chrom_lens = np.array([len(chroms[c]) for c in chrom_names])
    chrom_weights = chrom_lens / chrom_lens.sum()

    bio_seqs = []
    allowed = set("ACGT")
    while len(bio_seqs) < N_BIO:
        ci = rng.choice(len(chrom_names), p=chrom_weights)
        c = chrom_names[ci]
        pos = rng.integers(0, len(chroms[c]) - SEQ_LEN)
        s = chroms[c][pos:pos + SEQ_LEN]
        if len(s) == SEQ_LEN and all(ch in allowed for ch in s):
            bio_seqs.append(s)
    print(f"generated {len(bio_seqs)} hg38 sequences")

    all_seqs = np.array(dirichlet_seqs + bio_seqs)
    rng.shuffle(all_seqs)

    assert len(all_seqs) == N_TOTAL
    bad = sum(1 for s in all_seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    print(f"unique: {len(set(all_seqs))}")
    with open(OUT, "w") as fh:
        for s in all_seqs:
            fh.write(s + "\n")
    print(f"wrote {len(all_seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
