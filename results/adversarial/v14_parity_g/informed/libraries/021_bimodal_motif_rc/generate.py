#!/usr/bin/env python3
"""
021_bimodal_motif_rc — Combine 004 (bimodal motif vs random) + 019 (RC).
12.5k motif-packed + 12.5k random, each pair augmented with its RC. 50k total.
"""
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "sequences_0.txt"

L = 200
N_QUART = 12_500
SEED = 0

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
ACGT = "ACGT"
IUPAC = {"N": "ACGT", "S": "GC", "R": "AG", "Y": "CT", "W": "AT",
         "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
RC = str.maketrans("ACGT", "TGCA")


def resolve_iupac(s, rng):
    return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in s)


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


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
    motifs_clean = [resolve_iupac(m, rng) for m in MOTIFS]
    out = []
    for _ in range(N_QUART):
        s = motif_packed(rng, motifs_clean)
        out.append(s)
        out.append(s.translate(RC)[::-1])
    for _ in range(N_QUART):
        s = bg(rng, L)
        out.append(s)
        out.append(s.translate(RC)[::-1])
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} (motif+RC/random+RC) -> {OUT}")


if __name__ == "__main__":
    main()
