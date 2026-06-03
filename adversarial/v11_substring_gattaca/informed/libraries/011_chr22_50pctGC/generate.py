"""
E11: chr22_50pctGC

Sample 50k random 200bp windows from human chr22 (hg38), filtered to
~50% GC and no N's. Tests whether REAL biological sequences (with
natural dinuc/motif structure) help or hurt vs random_uniform.

Per-position distribution: matches natural human genome, which is
slightly AT-rich (~41% GC overall), with CpG depletion and other biases.
After filtering to 45-55% GC, dinuc structure should still differ from
uniform.

Prediction (per E1, E6 results): real biology will hurt by 0.05+ because
it's off-distribution for the IID-random eval. But a definitive test.
"""
import numpy as np
import os
import gzip

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
GC_LOW, GC_HIGH = 0.45, 0.55  # filter to match random_uniform's natural range
CHR22_PATH = os.path.join(os.path.dirname(__file__), "..", "..", "data", "chr22.fa.gz")

rng = np.random.default_rng(SEED)

def load_chr22():
    """Read chr22 FASTA, uppercase, concatenate."""
    with gzip.open(CHR22_PATH, "rt") as f:
        f.readline()  # header
        seq = "".join(line.strip().upper() for line in f)
    return seq

def main():
    print("Loading chr22...")
    seq = load_chr22()
    print(f"chr22 length: {len(seq):,}")

    out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
    sequences = []
    seq_len = len(seq)
    max_attempts = N_SEQS * 100
    attempts = 0
    while len(sequences) < N_SEQS and attempts < max_attempts:
        attempts += 1
        start = rng.integers(0, seq_len - SEQ_LEN + 1)
        window = seq[start : start + SEQ_LEN]
        if "N" in window:
            continue
        gc = (window.count("G") + window.count("C")) / SEQ_LEN
        if not (GC_LOW <= gc <= GC_HIGH):
            continue
        sequences.append(window)
    print(f"Collected {len(sequences)} sequences after {attempts} attempts")
    assert len(sequences) == N_SEQS

    with open(out_path, "w") as f:
        for s in sequences:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS} sequences to {out_path}")

if __name__ == "__main__":
    main()
