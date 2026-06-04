"""Experiment 015 — Motif-tiled sequences (structured non-random subspace).

Each 200bp sequence = several JASPAR motifs concatenated end-to-end, padded
to length 200 with small random spacers. Unlike 002 (3 motifs implanted on
180bp random scaffold ≈ 7-10% motif coverage), this packs the sequence
densely with motifs (~60-90% motif coverage).

Hypothesis (theory v8): the 0.397 plateau may be the ceiling for the GC-50
random subspace. A structured, motif-dense library lives in a fundamentally
different sequence subspace and may reach >0.402 if structure carries
transferable grammar — or fall sharply if structured sequences are
unrealistic for the model+evaluator.

Design: filter JASPAR motifs to length 6-15bp. For each sequence, sample
motifs until ~200bp of motif content reached, then fit them in random order,
with 1-3bp random spacers between motifs, padding to 200bp.

GC composition will be set by the motif distribution; we accept whatever
naturally emerges and report it.
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
from utils.seqlib import write_sequences, SEQLEN, is_clean

OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N = 50_000
MOTIF_MIN_LEN = 6
MOTIF_MAX_LEN = 15
SPACER_MIN = 1
SPACER_MAX = 3
BASES = "ACGT"


def random_spacer(rng, length):
    return "".join(rng.choice(list(BASES), size=length))


def make_tiled_seq(rng, motif_rng, pwms_pool):
    """Build a 200bp sequence by tiling motifs with short random spacers."""
    pieces = []
    cur_len = 0
    # Keep adding motif + spacer until we'd overflow 200, then pad.
    while cur_len < SEQLEN:
        # Sample motif from pool
        idx = int(motif_rng.integers(0, len(pwms_pool)))
        pwm = pwms_pool[idx]
        motif = sample_motif_instance(pwm, motif_rng)
        # Sample spacer
        sp_len = int(rng.integers(SPACER_MIN, SPACER_MAX + 1))
        # If adding this would overflow, decide:
        proj = cur_len + len(motif) + sp_len
        if proj > SEQLEN:
            # If motif alone still fits, add it
            if cur_len + len(motif) <= SEQLEN:
                pieces.append(motif)
                cur_len += len(motif)
            break
        pieces.append(motif)
        pieces.append(random_spacer(rng, sp_len))
        cur_len += len(motif) + sp_len
    seq = "".join(pieces)
    # Pad remaining
    if len(seq) < SEQLEN:
        seq += random_spacer(rng, SEQLEN - len(seq))
    seq = seq[:SEQLEN]
    return seq


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    motif_rng = np.random.default_rng(SEED + 7)

    print("Parsing JASPAR PFMs...")
    motifs = parse_jaspar()
    print(f"  Got {len(motifs)} PFMs")

    # Filter by length
    pwms_pool = [pwm for _, _, pwm in motifs
                 if MOTIF_MIN_LEN <= pwm.shape[1] <= MOTIF_MAX_LEN]
    print(f"  After length filter [{MOTIF_MIN_LEN}, {MOTIF_MAX_LEN}]: {len(pwms_pool)}")

    seqs = []
    for i in range(N):
        s = make_tiled_seq(rng, motif_rng, pwms_pool)
        assert is_clean(s), f"bad seq at {i}: {s[:30]}..."
        seqs.append(s)
        if (i + 1) % 10000 == 0:
            print(f"  generated {i+1}/{N}")

    # Stats
    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in seqs])
    print(f"GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}, min={gcs.min():.3f}, max={gcs.max():.3f}")

    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
