#!/usr/bin/env python3
"""
024_flagship_v2 — Flagship v1 minus Sharpr (009/013/014 showed Sharpr
provides no measurable advantage on v14), reallocated to more cCRE +
motif-packed sequences. RC pairs across the entire library.

Composition (50k total):
  - 12k cCRE PLS, 12k dELS, 4k pELS, 4k CTCF-only, 4k DNase-H3K4me3
  - 7k synthetic motif-packed (50-TF vocab)
  - 7k RC of those motif-packed
Total = 36k cCRE + 14k motif(+RC).
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
SEED = 0
ACGT_SET = set("ACGT")
ACGT = "ACGT"

PARAM = {"PLS": 12_000, "dELS": 12_000, "pELS": 4_000,
         "CTCF-only": 4_000, "DNase-H3K4me3": 4_000}
N_MOTIF = 7_000
N_MOTIF_RC = 7_000

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
RC = str.maketrans("ACGT", "TGCA")


def resolve_iupac(s, rng):
    return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in s)


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


def base_class(label):
    for p in label.split(","):
        if p in PARAM:
            return p
    return None


def motif_packed(rng, motifs_clean):
    n_motifs = rng.randint(4, 8)
    seq = list(bg(rng, L))
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
    return "".join(seq)


def main():
    rng = random.Random(SEED)
    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)
    motifs_clean = [resolve_iupac(m, rng) for m in MOTIFS]

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

    out = []
    for cls, n in PARAM.items():
        taken = 0
        for chrom, mid in by_class[cls]:
            if taken >= n:
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
            out.append(seq)
            taken += 1
        print(f"  {cls}: {taken}")

    motifs_seqs = [motif_packed(rng, motifs_clean) for _ in range(N_MOTIF)]
    out += motifs_seqs
    out += [s.translate(RC)[::-1] for s in motifs_seqs[:N_MOTIF_RC]]
    print(f"  motif: {N_MOTIF}, motif-RC: {N_MOTIF_RC}")

    assert len(out) == 50_000, f"got {len(out)}"
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} -> {OUT}")


if __name__ == "__main__":
    main()
