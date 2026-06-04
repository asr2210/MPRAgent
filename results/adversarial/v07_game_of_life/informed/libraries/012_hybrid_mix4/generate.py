"""Experiment 012 — Mixed-strategy hybrid library (4-way blend).

50k sequences = 12.5k from each of 4 known-good single strategies:
1. random + Malinois top (12.5k from 007's pool)
2. random + Malinois span (12.5k from 006's pool)
3. cCRE GC-matched (12.5k fresh sample)
4. gc_50 i.i.d. random (12.5k fresh)

All sub-pools have GC ≈ 0.50.

Hypothesis: combining multiple complementary distributions hedges against
unknown eval-set bias. If the trained model benefits from diverse forms of
training signal, this should exceed the best single strategy (007 = 0.397).
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
from utils.seqlib import extract, write_sequences, fasta, SEQLEN

CCRE_BED = os.path.join(REPO, "data", "ENCODE_cCREs_v3.bed")
OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
N_TARGET = 50_000
PER = N_TARGET // 4   # 12500
GC_MIN, GC_MAX = 0.45, 0.55


def load_sequences_from_file(path, n, rng):
    """Load first n lines from a sequence file; shuffle within if requested."""
    with open(path) as f:
        seqs = [line.strip() for line in f if line.strip()]
    rng.shuffle(seqs)
    return seqs[:n]


def gen_gc50_random(rng, n):
    """Generate n random GC-50 200bp ACGT sequences."""
    lookup = np.array(list("ACGT"))
    idx = rng.integers(0, 4, size=(n, SEQLEN), dtype=np.int8)
    return ["".join(row) for row in lookup[idx]]


def gen_ccre_gc_matched(rng, n):
    """Sample cCRE-derived 200bp windows with GC in [0.45, 0.55]."""
    fa = fasta()
    print(f"  Loading cCREs from {CCRE_BED}...")
    ccres = []
    with open(CCRE_BED) as f:
        for line in f:
            p = line.rstrip().split("\t")
            if len(p) < 3:
                continue
            try:
                ccres.append((p[0], int(p[1]), int(p[2])))
            except ValueError:
                continue
    rng.shuffle(ccres)
    print(f"  Have {len(ccres):,} cCREs")
    seqs = []
    seen = set()
    for chrom, s, e in ccres:
        if len(seqs) >= n:
            break
        center = (s + e) // 2
        for off in (0, -50, 50, -100, 100):
            seq = extract(chrom, center + off)
            if seq is None or seq in seen:
                continue
            gc = (seq.count("G") + seq.count("C")) / SEQLEN
            if gc < GC_MIN or gc > GC_MAX:
                continue
            seen.add(seq)
            seqs.append(seq)
            break
    return seqs[:n]


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)

    print(f"Sampling {PER:,} from 007 random+Malinois top...")
    p007 = os.path.join(REPO, "libraries", "007_malinois_top_magnitude", "sequences_0.txt")
    s007 = load_sequences_from_file(p007, PER, rng)
    print(f"  got {len(s007)}")

    print(f"Sampling {PER:,} from 006 random+Malinois span...")
    p006 = os.path.join(REPO, "libraries", "006_malinois_oracle_activity_span", "sequences_0.txt")
    s006 = load_sequences_from_file(p006, PER, rng)
    print(f"  got {len(s006)}")

    print(f"Generating {PER:,} fresh cCRE GC-matched...")
    s_ccre = gen_ccre_gc_matched(rng, PER)
    print(f"  got {len(s_ccre)}")

    print(f"Generating {PER:,} fresh GC-50 random...")
    s_rand = gen_gc50_random(rng, PER)
    print(f"  got {len(s_rand)}")

    all_seqs = s007 + s006 + s_ccre + s_rand
    # De-duplicate (preserving randomness)
    seen = set()
    uniq = []
    for s in all_seqs:
        if s not in seen and len(s) == SEQLEN:
            seen.add(s)
            uniq.append(s)
    print(f"  Total after dedup: {len(uniq):,}")
    if len(uniq) < N_TARGET:
        # Top up with fresh random
        extra = gen_gc50_random(rng, N_TARGET - len(uniq) + 100)
        for s in extra:
            if s not in seen:
                seen.add(s)
                uniq.append(s)
                if len(uniq) >= N_TARGET:
                    break
    rng.shuffle(uniq)
    uniq = uniq[:N_TARGET]

    gcs = np.array([(s.count("G") + s.count("C")) / SEQLEN for s in uniq])
    print(f"Final GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")
    write_sequences(uniq, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
