#!/usr/bin/env python3
"""
004_bimodal_motif_vs_random — Bimodal library.
- 25,000 sequences are motif-packed (4-8 strong TF motifs each).
- 25,000 sequences are random uniform 200bp (no motif insertion).
Interleaved to mix during training.

Tests whether the v14 model can learn a "motif present → high activity"
contrast when given strong bimodal label variance. If yes, r should jump
above the noise floor.
"""
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "sequences_0.txt"

L = 200
N_HALF = 25_000
SEED = 0

MOTIFS = [
    "TGACTCA", "TGACGTCA", "GGGACTTTCC", "GGGGCGGGG", "TATAAA",
    "CACGTG", "CCAATCA", "GATAAG", "AGGAAG", "TTGCGCAAT",
    "TGAACTTTGAACT", "GTTAATCATTAAC", "CAGCTG", "GGGGAGGGG",
    "AACAAAG", "AGGTCA", "TGTGGTTT",
]

ACGT = "ACGT"


def bg(rng, n):
    return "".join(rng.choices(ACGT, k=n))


def motif_packed(rng):
    n_motifs = rng.randint(4, 8)
    placements = []
    occupied = []
    for _ in range(n_motifs):
        m = rng.choice(MOTIFS)
        for _ in range(20):
            start = rng.randint(0, L - len(m))
            end = start + len(m)
            if all(end <= s or start >= e for s, e in occupied):
                placements.append((start, m))
                occupied.append((start, end))
                break
    seq = list(bg(rng, L))
    for start, m in placements:
        for i, ch in enumerate(m):
            seq[start + i] = ch
    return "".join(seq)


def random_seq(rng):
    return bg(rng, L)


def main():
    rng = random.Random(SEED)
    motif_seqs = [motif_packed(rng) for _ in range(N_HALF)]
    rand_seqs = [random_seq(rng) for _ in range(N_HALF)]
    # Interleave
    all_seq = []
    for a, b in zip(motif_seqs, rand_seqs):
        all_seq.append(a)
        all_seq.append(b)
    assert len(all_seq) == 50_000
    with open(OUT, "w") as f:
        for s in all_seq:
            f.write(s + "\n")
    print(f"Wrote 25k motif + 25k random -> {OUT}")


if __name__ == "__main__":
    main()
