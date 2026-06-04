"""
Experiment 002: uniform DHS sampling.

Hypothesis: pure uniform random sampling over 3.59M DHS elements should
reproduce dhs_random baseline (~0.7089 eval_01). This calibrates my
pipeline against the prior numbers.

If the result is far off ~0.71, the issue is in my pipeline (FASTA, summit,
window, ATCG handling). If close, my pipeline is sound and the prior
mean_signal weighting was the problem.
"""
import gzip
import os
import numpy as np
import twobitreader

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DHS_FILE = os.path.join(ROOT, "data", "DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz")
TWOBIT = os.path.join(ROOT, "data", "hg38.2bit")
OUT = os.path.join(HERE, "sequences_0.txt")

WINDOW = 200
N_SEQ = 50_000
SEED = 1
CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}


def load_dhs():
    rows = []
    with gzip.open(DHS_FILE, "rt") as f:
        f.readline()  # header
        for line in f:
            p = line.rstrip("\n").split("\t")
            if p[0] not in CHROMS:
                continue
            rows.append((p[0], int(p[6])))  # chrom, summit
    return rows


def extract(rows, idxs, genome):
    half = WINDOW // 2
    out = []
    for i in idxs:
        chrom, summit = rows[i]
        start = max(0, summit - half)
        end = start + WINDOW
        seq = genome[chrom][start:end].upper()
        if len(seq) == WINDOW and all(c in "ACGT" for c in seq):
            out.append(seq)
    return out


def main():
    rows = load_dhs()
    print(f"DHS rows: {len(rows):,}")

    rng = np.random.default_rng(SEED)
    # Oversample to allow rejecting N-containing windows.
    idxs = rng.choice(len(rows), size=int(N_SEQ * 1.3), replace=False)
    genome = twobitreader.TwoBitFile(TWOBIT)
    seqs = extract(rows, idxs, genome)
    print(f"After N-filter: {len(seqs)}")

    if len(seqs) < N_SEQ:
        extra = rng.choice(len(rows), size=int(N_SEQ * 0.5), replace=False)
        seqs += extract(rows, extra, genome)
    seqs = seqs[:N_SEQ]

    assert len(seqs) == N_SEQ
    with open(OUT, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQ} sequences")


if __name__ == "__main__":
    main()
