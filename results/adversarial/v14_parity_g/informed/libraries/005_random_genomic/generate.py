#!/usr/bin/env python3
"""
005_random_genomic — 50,000 random 200bp windows uniformly drawn from across
the hg38 main assembly (not curated to cCREs). Compare to 001 (cCRE-curated).

Test: does broad genomic sampling (mostly intergenic / inactive) work better
than curated regulatory regions? If v14 cares about overall genome composition
rather than enhancer biology, this should beat 001. If not, biology helps but
we still might be in the noise floor.
"""
import random
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data"
GENOME = DATA / "hg38.fa"
OUT = HERE / "sequences_0.txt"

N_SEQ = 50_000
L = 200
SEED = 0
ACGT = set("ACGT")

# Use only the 22 autosomes + X + Y (skip chrM and unplaced contigs).
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def main():
    rng = random.Random(SEED)
    np.random.seed(SEED)
    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)

    chrom_lens = {c: len(fa[c]) for c in CHROMS if c in fa.keys()}
    total = sum(chrom_lens.values())
    chrom_p = [chrom_lens[c] / total for c in chrom_lens]
    chrom_list = list(chrom_lens.keys())

    sequences = []
    attempts = 0
    while len(sequences) < N_SEQ and attempts < N_SEQ * 50:
        attempts += 1
        c = np.random.choice(chrom_list, p=chrom_p)
        clen = chrom_lens[c]
        s = rng.randint(0, clen - L - 1)
        seq = str(fa[c][s:s + L]).upper()
        if len(seq) != L or not set(seq).issubset(ACGT):
            continue
        sequences.append(seq)
    print(f"Collected {len(sequences)} sequences in {attempts} attempts")
    assert len(sequences) == N_SEQ
    with open(OUT, "w") as f:
        for s in sequences:
            f.write(s + "\n")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
