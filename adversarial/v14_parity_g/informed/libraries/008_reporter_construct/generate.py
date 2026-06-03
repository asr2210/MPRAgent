#!/usr/bin/env python3
"""
008_reporter_construct — Sequences that mimic the structure of a typical
synthetic MPRA reporter construct: enhancer cassette + minimal promoter.

Structure (200bp):
  positions   0 - 140  : variable enhancer region with 4-8 TF motifs
  positions 140 - 200  : fixed 60bp minimal-promoter-like cassette
                         (TATA box + Inr-like + downstream spacer)

Test: does giving the v14 model a stereotyped construct structure (the kind
of structure most MPRA training data has) let it learn anything?
"""
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "sequences_0.txt"

L = 200
N_SEQ = 50_000
SEED = 0

# Fixed minP-like cassette (60bp).
# Approximation: GC-balanced spacer + TATA + Inr + downstream spacer.
MINP = "GCAGCAGCGGTTTCGCTAGCTATAAAAGGCCGTCAGTACTCAGCTGAACGTCAGCTCAAA"
assert len(MINP) == 60

ENH_LEN = L - len(MINP)  # 140
assert ENH_LEN == 140

MOTIFS = [
    "TGACTCA", "TGACGTCA", "GGGACTTTCC", "GGGGCGGGG",
    "CACGTG", "CCAATCA", "GATAAG", "AGGAAG",
    "TTGCGCAAT", "CAGCTG", "GGGGAGGGG", "AGGTCA",
    "TGTGGTTT", "AACAAAG",
]
ACGT = "ACGT"


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


def enhancer_with_motifs(rng):
    n_motifs = rng.randint(4, 8)
    occupied = []
    placements = []
    for _ in range(n_motifs):
        m = rng.choice(MOTIFS)
        for _ in range(20):
            start = rng.randint(0, ENH_LEN - len(m))
            end = start + len(m)
            if all(end <= s or start >= e for s, e in occupied):
                placements.append((start, m))
                occupied.append((start, end))
                break
    seq = list(bg(rng, ENH_LEN))
    for start, m in placements:
        for i, ch in enumerate(m):
            seq[start + i] = ch
    return "".join(seq)


def main():
    rng = random.Random(SEED)
    with open(OUT, "w") as f:
        for _ in range(N_SEQ):
            enh = enhancer_with_motifs(rng)
            f.write(enh + MINP + "\n")
    print(f"Wrote {N_SEQ} reporter-construct sequences -> {OUT}")


if __name__ == "__main__":
    main()
