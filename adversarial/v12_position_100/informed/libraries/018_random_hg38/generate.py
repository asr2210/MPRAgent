#!/usr/bin/env python3
"""
Experiment 018 — Random hg38 200bp regions.

Sample 50k random 200bp windows from hg38 autosomes + chrX. Reject windows
containing N. Tests whether generic genomic (non-DHS) sequences give different
eval performance than DHS-enriched biology.
"""
import os
import gzip
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
HG38 = os.path.abspath(os.path.join(HERE, "..", "..", "data", "hg38.fa.gz"))
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 0
N_TOTAL = 50_000
SEQ_LEN = 200
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
                    print(f"loaded {cur_name}: {len(chroms[cur_name]):,} bp")
                cur_name = line[1:].split()[0]
                cur_chunks = []
                if cur_name not in USE_CHROMS and len(chroms) == len(USE_CHROMS):
                    break
            else:
                if cur_name in USE_CHROMS:
                    cur_chunks.append(line)
        if cur_name in USE_CHROMS and cur_name not in chroms:
            chroms[cur_name] = "".join(cur_chunks).upper()
            print(f"loaded {cur_name}: {len(chroms[cur_name]):,} bp")
    return chroms


def main():
    rng = np.random.default_rng(SEED)
    chroms = load_chroms()
    chrom_names = sorted(chroms.keys())
    chrom_lens = np.array([len(chroms[c]) for c in chrom_names])
    chrom_weights = chrom_lens / chrom_lens.sum()
    print(f"total: {chrom_lens.sum():,} bp across {len(chrom_names)} chroms")

    seqs = []
    allowed = set("ACGT")
    attempts = 0
    while len(seqs) < N_TOTAL:
        # Sample chromosome weighted by length
        n_batch = min(N_TOTAL * 2, (N_TOTAL - len(seqs)) * 3)
        chrom_idx = rng.choice(len(chrom_names), size=n_batch, p=chrom_weights)
        for ci in chrom_idx:
            if len(seqs) >= N_TOTAL:
                break
            c = chrom_names[ci]
            pos = rng.integers(0, len(chroms[c]) - SEQ_LEN)
            s = chroms[c][pos:pos + SEQ_LEN]
            if len(s) == SEQ_LEN and all(ch in allowed for ch in s):
                seqs.append(s)
            attempts += 1
        if attempts > N_TOTAL * 20:
            raise RuntimeError(f"too many rejections; got {len(seqs)}/{N_TOTAL}")

    seqs = np.array(seqs)
    rng.shuffle(seqs)
    print(f"sampled {len(seqs)} clean sequences from {attempts} attempts")

    assert len(seqs) == N_TOTAL
    bad = sum(1 for s in seqs if len(s) != SEQ_LEN or any(c not in allowed for c in s))
    assert bad == 0, f"{bad} bad sequences"
    print(f"unique: {len(set(seqs))}")
    with open(OUT, "w") as fh:
        for s in seqs:
            fh.write(s + "\n")
    print(f"wrote {len(seqs)} sequences to {OUT}")


if __name__ == "__main__":
    main()
