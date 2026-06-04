#!/usr/bin/env python3
"""
Experiment 026 — Chaotic 5-way mix of best strategies.

10k each from: Dirichlet(0.3), Dirichlet(0.1), shuffled DHS, hg38 random,
Markov-Dirichlet(0.3). Tests if maximum strategy diversity helps.
"""
import os
import gzip
import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
DHS_DATA = os.path.abspath(os.path.join(HERE, "..", "..", "data", "train_all_classifier_light.csv.gz"))
HG38 = os.path.abspath(os.path.join(HERE, "..", "..", "data", "hg38.fa.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_PER = 10_000
SEQ_LEN = 200
ALPHABET = np.array(list("ACGT"))
USE_CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX"}


def gen_dirichlet(rng, n, alpha):
    seqs = []
    for _ in range(n):
        p = rng.dirichlet([alpha] * 4)
        idx = rng.choice(4, size=SEQ_LEN, p=p)
        seqs.append("".join(ALPHABET[idx]))
    return seqs


def gen_shuffled_dhs(rng, n):
    df = pd.read_csv(DHS_DATA, sep="\t")
    pool = df["raw_sequence"].values
    idx = rng.choice(len(pool), size=n, replace=False)
    seqs = []
    for s in pool[idx]:
        arr = np.array(list(s))
        rng.shuffle(arr)
        seqs.append("".join(arr))
    return seqs


def gen_hg38(rng, n):
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
    chrom_names = sorted(chroms.keys())
    chrom_lens = np.array([len(chroms[c]) for c in chrom_names])
    chrom_weights = chrom_lens / chrom_lens.sum()
    seqs = []
    allowed = set("ACGT")
    while len(seqs) < n:
        ci = rng.choice(len(chrom_names), p=chrom_weights)
        c = chrom_names[ci]
        pos = rng.integers(0, len(chroms[c]) - SEQ_LEN)
        s = chroms[c][pos:pos + SEQ_LEN]
        if len(s) == SEQ_LEN and all(ch in allowed for ch in s):
            seqs.append(s)
    return seqs


def gen_markov_dirichlet(rng, n):
    seqs = []
    for _ in range(n):
        T = np.array([rng.dirichlet([0.3] * 4) for _ in range(4)])
        # start uniform
        cur = rng.integers(4)
        chars = [ALPHABET[cur]]
        for _ in range(SEQ_LEN - 1):
            cur = rng.choice(4, p=T[cur])
            chars.append(ALPHABET[cur])
        seqs.append("".join(chars))
    return seqs


def main():
    rng = np.random.default_rng(SEED)
    parts = []
    print("gen Dirichlet(0.3)...")
    parts.extend(gen_dirichlet(rng, N_PER, 0.3))
    print("gen Dirichlet(0.1)...")
    parts.extend(gen_dirichlet(rng, N_PER, 0.1))
    print("gen shuffled DHS...")
    parts.extend(gen_shuffled_dhs(rng, N_PER))
    print("gen hg38...")
    parts.extend(gen_hg38(rng, N_PER))
    print("gen Markov-Dirichlet...")
    parts.extend(gen_markov_dirichlet(rng, N_PER))

    seqs = np.array(parts)
    rng.shuffle(seqs)
    assert len(seqs) == 50_000
    allowed = set("ACGT")
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    print(f"unique: {len(set(seqs))}")
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
