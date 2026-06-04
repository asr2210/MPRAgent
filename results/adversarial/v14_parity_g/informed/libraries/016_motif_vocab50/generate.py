#!/usr/bin/env python3
"""
016_motif_vocab50 — Same motif-packed strategy as 003, but with a 50-TF
vocabulary (vs 14 in 003/012). Each sequence carries 4-8 motifs from a
broader TF panel covering more regulatory programs.

Rationale: if the model is motif-aware, a richer motif vocabulary expands
the regulatory grammar it can learn. Reasonable lift over the 14-motif test.
"""
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "sequences_0.txt"

L = 200
N_SEQ = 50_000
SEED = 0

# 50 TF consensus motifs spanning AP-1, CREB, NF-kB, SP/KLF, E-box, GATA,
# ETS, C/EBP, HNF, Fox, Sox, Oct, IRF, STAT, RUNX, TEAD, POU, SMAD, p53,
# nuclear receptors, etc.
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


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


def resolve_iupac(s, rng):
    # Resolve N/S/R/Y/etc. ambiguity codes to random ACGT
    iupac = {"N": "ACGT", "S": "GC", "R": "AG", "Y": "CT", "W": "AT",
             "K": "GT", "M": "AC", "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG"}
    return "".join(rng.choice(iupac[c]) if c in iupac else c for c in s)


def main():
    rng = random.Random(SEED)
    motifs_clean = [resolve_iupac(m, rng) for m in MOTIFS]  # static resolution

    with open(OUT, "w") as f:
        for _ in range(N_SEQ):
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
            f.write("".join(seq) + "\n")
    print(f"Wrote {N_SEQ} 50-vocab motif-packed sequences -> {OUT}")


if __name__ == "__main__":
    main()
