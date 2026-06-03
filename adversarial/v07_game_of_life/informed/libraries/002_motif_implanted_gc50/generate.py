"""Experiment 002 — GC-50 random scaffold with implanted JASPAR motifs.

Hypothesis: if biological elements don't help (exp 001), but GC=50% random does
(strategies.md baseline ~0.40), perhaps the model can learn motif content if
delivered on a clean, composition-matched scaffold. Implant 3 random JASPAR
motifs per 200bp GC-50 random sequence at non-overlapping random positions.

Expected score: > 0.40 if motif content helps the model; ~0.40 if it doesn't.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
from utils.jaspar import parse_jaspar, sample_motif_instance
from utils.seqlib import write_sequences, is_clean

SEED = 1
N = 50000
SEQLEN = 200
MOTIFS_PER_SEQ = 3
GC_TARGET = 0.5


def random_gc50(rng, length=SEQLEN):
    """Generate i.i.d. random sequence with each base equally likely (gives ~50% GC)."""
    return rng.choice(list("ACGT"), size=length)


def implant_motifs(scaffold_arr, motif_strs, rng):
    """Implant motifs at random non-overlapping positions in the scaffold (numpy array)."""
    L = len(scaffold_arr)
    placed = []  # list of (start, end)
    out = scaffold_arr.copy()
    tries = 0
    for motif in motif_strs:
        ml = len(motif)
        if ml >= L:
            continue
        for _ in range(30):
            start = int(rng.integers(0, L - ml + 1))
            end = start + ml
            if all(end <= s or start >= e for s, e in placed):
                placed.append((start, end))
                out[start:end] = list(motif)
                break
        tries += 1
    return out


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    motifs = parse_jaspar()
    print(f"Parsed {len(motifs)} JASPAR PFMs")
    pwms = [pwm for _, _, pwm in motifs]

    seqs = []
    motif_rng = np.random.default_rng(SEED + 7)
    for i in range(N):
        scaffold = random_gc50(rng)
        # Pick MOTIFS_PER_SEQ random motifs and sample instances
        motif_indices = motif_rng.choice(len(pwms), size=MOTIFS_PER_SEQ, replace=False)
        motif_strs = [sample_motif_instance(pwms[j], motif_rng) for j in motif_indices]
        out_arr = implant_motifs(scaffold, motif_strs, motif_rng)
        seq = "".join(out_arr)
        assert is_clean(seq), f"bad seq at {i}: {seq[:20]}..."
        seqs.append(seq)
        if (i + 1) % 10000 == 0:
            print(f"  generated {i+1}/{N}")

    out_path = os.path.join(THIS, "sequences_0.txt")
    write_sequences(seqs, out_path)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
