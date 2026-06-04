#!/usr/bin/env python3
"""
018_vocab50_dense — 016 design (50-TF motif vocab) but denser: 8-16 motifs
per sequence instead of 4-8. Tests whether more motif occupancy strengthens
signal.
"""
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "sequences_0.txt"

L = 200
N_SEQ = 50_000
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


def resolve_iupac(s, rng):
    return "".join(rng.choice(IUPAC[c]) if c in IUPAC else c for c in s)


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


def main():
    rng = random.Random(SEED)
    motifs_clean = [resolve_iupac(m, rng) for m in MOTIFS]
    with open(OUT, "w") as f:
        for _ in range(N_SEQ):
            n_motifs = rng.randint(8, 16)
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
            f.write("".join(seq) + "\n")
    print(f"Wrote {N_SEQ} dense motif-packed sequences -> {OUT}")


if __name__ == "__main__":
    main()
