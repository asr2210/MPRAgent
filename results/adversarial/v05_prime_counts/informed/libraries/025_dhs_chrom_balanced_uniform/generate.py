"""
025_dhs_chrom_balanced_uniform
==============================
33/67 mix (best ratio from Exps 016-022) but DHS half is CHROM-BALANCED:
~700 dhs_specific from each of 24 chromosomes (~16.7k bio) + 33.3k uniform.

Tests if balanced chrom representation in the bio half pushes higher than
weighted random sampling (which preferentially hits big chroms).
"""
import gzip
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_BIO = 16_667
N_UNIFORM = 50_000 - N_BIO
N_SEQS = 50_000
SEQ_LEN = 200
HALF = SEQ_LEN // 2
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
DHS_PATH = DATA_DIR / "dhs_index.txt.gz"
HG38_PATH = DATA_DIR / "hg38.fa"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
PER_CHROM = N_BIO // len(CHROMS)  # 694
EXTRA = N_BIO - PER_CHROM * len(CHROMS)


def load_dhs():
    by_chrom = defaultdict(list)
    with gzip.open(DHS_PATH, "rt") as f:
        f.readline()
        for line in f:
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            if chrom not in CHROMS:
                continue
            by_chrom[chrom].append((int(parts[6]), int(parts[5])))
    return by_chrom


def main():
    rng = np.random.default_rng(SEED)
    by_chrom = load_dhs()
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)

    bio_seqs = []
    for ci, chrom in enumerate(CHROMS):
        target = PER_CHROM + (1 if ci < EXTRA else 0)
        sites = by_chrom[chrom]
        nsamp = np.array([s[1] for s in sites])
        w = 1.0 / np.sqrt(nsamp.astype(np.float64))
        w /= w.sum()
        pool = rng.choice(len(sites), size=int(target * 1.5), replace=False, p=w)
        chr_seqs = []
        for i in pool:
            if len(chr_seqs) >= target:
                break
            summit = sites[i][0]
            L = len(fa[chrom])
            start = max(0, summit - HALF)
            end = start + SEQ_LEN
            if end > L:
                end = L
                start = end - SEQ_LEN
            if start < 0:
                continue
            s = str(fa[chrom][start:end])
            if len(s) != SEQ_LEN or set(s) - set("ACGT"):
                continue
            chr_seqs.append(s)
        bio_seqs.extend(chr_seqs)
        print(f"  {chrom}: {len(chr_seqs)}/{target}", file=sys.stderr)
    assert len(bio_seqs) == N_BIO, len(bio_seqs)

    bases = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(N_UNIFORM, SEQ_LEN))
    uniform_seqs = ["".join(row) for row in bases[idx]]

    all_seqs = bio_seqs + uniform_seqs
    rng.shuffle(all_seqs)
    out_path = Path(__file__).parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        for s in all_seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
