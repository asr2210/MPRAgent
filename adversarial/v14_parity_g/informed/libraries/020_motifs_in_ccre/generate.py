#!/usr/bin/env python3
"""
020_motifs_in_ccre — Use real ENCODE cCRE genomic sequences as backbone,
then insert 4-8 TF motifs from the 50-vocab into each. Combines real
genomic local statistics with synthetic combinatorial motif coverage.

Tests: (real genomic backbone) + (controlled motif insertion) might
provide both authentic flanking-base statistics AND controlled motif
combinations, which neither 005/006 (real-only) nor 016 (random+motifs)
captures.
"""
import random
from collections import defaultdict
from pathlib import Path

from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
CCRE = ROOT / "data" / "GRCh38-cCREs.bed"
GENOME = ROOT / "data" / "hg38.fa"
OUT = HERE / "sequences_0.txt"

L = 200
N_SEQ = 50_000
SEED = 0
ACGT_SET = set("ACGT")
ACGT = "ACGT"

MOTIFS = [
    "TGACTCA", "TGAGTCA", "TGACGTCA", "GGGACTTTCC", "GGGAATTCCC",
    "GGGGCGGGG", "GGGGGCGGGG", "CACGTG", "CAGGTG", "CAGCTG",
    "CACGCG", "CCAATCA", "ATTGGC", "GATAAG", "AGATAA",
    "GGAAGT", "AGGAAG", "CGGAAG", "TTGCGCAAT", "ATTGCGCAA",
    "GTTAATNATTAAC", "TTAATGA", "TGTGGTTT", "AACAAAG", "AGGTCA",
    "TGACCT", "AGAACA", "TGTTCT", "GTAAACA", "TGTTTAC",
    "AACAATG", "CATTGTT", "TAATTA", "TAATGG", "CCATTA",
    "GGGGAGGGG", "CCCCGCCC", "ATTTGCAT", "ATGCAAAT", "AANAGTGT",
    "ACACTTNNT", "TGASTCAGCA", "TTCNNNGAA", "TTCCGGGAA", "CAGGAAG",
    "ACCGGAAG", "TGTGGAAA", "TTTCCACA", "CATATG", "GCCNNNGGC",
]
IUPAC = {"N": "ACGT", "S": "GC", "R": "AG", "Y": "CT", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}


def resolve_iupac(s, rng):
    return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in s)


def base_class(label):
    for p in label.split(","):
        if p in ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"):
            return p
    return None


def main():
    rng = random.Random(SEED)
    motifs_clean = [resolve_iupac(m, rng) for m in MOTIFS]

    # Load balanced cCRE pool: 10k from each of 5 classes
    by_class = defaultdict(list)
    with open(CCRE) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            if "_" in chrom or chrom == "chrM":
                continue
            cls = base_class(parts[-1])
            if cls is None:
                continue
            try:
                start = int(parts[1]); end = int(parts[2])
            except ValueError:
                continue
            by_class[cls].append((chrom, (start + end) // 2))
    for cls in by_class:
        rng.shuffle(by_class[cls])

    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)
    backbones = []
    per_cls = N_SEQ // 5
    for cls in ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"):
        taken = 0
        for chrom, mid in by_class[cls]:
            if taken >= per_cls:
                break
            s = mid - L // 2
            e = s + L
            if s < 0:
                continue
            try:
                clen = len(fa[chrom])
            except KeyError:
                continue
            if e > clen:
                continue
            seq = str(fa[chrom][s:e]).upper()
            if len(seq) != L or not set(seq).issubset(ACGT_SET):
                continue
            backbones.append(seq)
            taken += 1
        print(f"  backbone {cls}: {taken}")

    out = []
    for bb in backbones:
        seq = list(bb)
        n_motifs = rng.randint(4, 8)
        occupied = []
        for _ in range(n_motifs):
            m = rng.choice(motifs_clean)
            for _ in range(20):
                start = rng.randint(0, L - len(m))
                end = start + len(m)
                if all(end <= s or start >= e for s, e in occupied):
                    for i, ch in enumerate(m):
                        seq[start + i] = ch
                    occupied.append((start, end))
                    break
        out.append("".join(seq))

    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} motif-in-cCRE sequences -> {OUT}")


if __name__ == "__main__":
    main()
