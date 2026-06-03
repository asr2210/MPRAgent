"""
006_random_genomic
==================
50k 200bp windows from random positions anywhere in hg38 (no annotation
filter). Mostly intergenic / repetitive / non-regulatory regions.

Hypothesis: if the eval favors unannotated genomic sequences (~99% of
genome), this should score higher than cCRE/DHS. If eval favors regulatory
elements specifically, this should score lower. Reveals whether eval cares
about regulatory content.
"""
import sys
from pathlib import Path

import numpy as np
from pyfaidx import Fasta

SEED = 0
N_SEQS = 50_000
SEQ_LEN = 200
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
HG38_PATH = DATA_DIR / "hg38.fa"
CANONICAL = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def main():
    rng = np.random.default_rng(SEED)
    fa = Fasta(str(HG38_PATH), as_raw=False, sequence_always_upper=True)
    chrom_lens = {c: len(fa[c]) for c in CANONICAL}
    total = sum(chrom_lens.values())
    # Sample chromosomes proportional to length, then a random position on each.
    chr_names = list(chrom_lens.keys())
    chr_probs = np.array([chrom_lens[c] for c in chr_names], dtype=np.float64)
    chr_probs /= chr_probs.sum()

    out_path = Path(__file__).parent / "sequences_0.txt"
    seqs = []
    while len(seqs) < N_SEQS:
        # Batch sample
        batch = max(1000, N_SEQS - len(seqs))
        chr_idx = rng.choice(len(chr_names), size=batch, p=chr_probs)
        for i in chr_idx:
            chrom = chr_names[i]
            L = chrom_lens[chrom]
            start = int(rng.integers(0, L - SEQ_LEN + 1))
            s = str(fa[chrom][start:start + SEQ_LEN])
            if len(s) != SEQ_LEN or set(s) - set("ACGT"):
                continue
            seqs.append(s)
            if len(seqs) >= N_SEQS:
                break
        if len(seqs) % 10_000 < 1000:
            print(f"  {len(seqs)}/{N_SEQS}", file=sys.stderr)

    with open(out_path, "w") as f:
        for s in seqs:
            f.write(s + "\n")
    print(f"Wrote {N_SEQS}", file=sys.stderr)


if __name__ == "__main__":
    main()
