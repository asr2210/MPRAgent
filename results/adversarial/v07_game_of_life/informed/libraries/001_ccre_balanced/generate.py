"""Experiment 001 — ENCODE cCRE category-balanced library.

Sample 10,000 sequences from each of 5 cCRE categories (50,000 total).
Categories: PLS, pELS, dELS, CTCF-only, DNase-H3K4me3 (CTCF-bound variants collapsed).
Each sequence is a 200bp window centered on the cCRE midpoint, extracted from hg38.

Hypothesis: balanced functional categories teach universal regulatory primitives that
generalize across cell types, beating cell-type-biased DHS sampling (dhs_topic 0.7232).
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
from utils.seqlib import extract_from_interval, write_sequences

CCRE_BED = os.path.join(REPO, "data", "ENCODE_cCREs_v3.bed")
OUT = os.path.join(THIS, "sequences_0.txt")
SEED = 1
PER_CAT = 10000  # 5 categories * 10k = 50k

CATEGORIES = {
    "PLS": ["PLS", "PLS,CTCF-bound"],
    "pELS": ["pELS", "pELS,CTCF-bound"],
    "dELS": ["dELS", "dELS,CTCF-bound"],
    "CTCF": ["CTCF-only,CTCF-bound"],
    "DNase-H3K4me3": ["DNase-H3K4me3", "DNase-H3K4me3,CTCF-bound"],
}


def load_ccres_by_category():
    """Return dict: category_label -> list of (chrom, start, end)."""
    by_cat = {k: [] for k in CATEGORIES}
    label_to_cat = {}
    for cat, labels in CATEGORIES.items():
        for lbl in labels:
            label_to_cat[lbl] = cat
    with open(CCRE_BED) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            if len(parts) < 6:
                continue
            chrom, start, end, _, _, label = parts[:6]
            cat = label_to_cat.get(label)
            if cat is None:
                continue
            by_cat[cat].append((chrom, int(start), int(end)))
    return by_cat


def sample_from_category(intervals, n, rng):
    """Randomly sample n valid 200bp sequences from a list of intervals.
    Each interval may produce one sequence; if extraction fails (N bases,
    boundary), draw another. Sample without replacement until we have n."""
    seqs = []
    seen = set()
    idxs = rng.permutation(len(intervals))
    pos = 0
    while len(seqs) < n and pos < len(idxs):
        i = idxs[pos]
        pos += 1
        chrom, s, e = intervals[i]
        seq = extract_from_interval(chrom, s, e)
        if seq is None or seq in seen:
            continue
        seen.add(seq)
        seqs.append(seq)
    if len(seqs) < n:
        # fall back: try jittered windows from already-tried intervals
        rng2 = np.random.default_rng(rng.integers(1 << 31))
        for _ in range(10):
            if len(seqs) >= n:
                break
            for i in rng2.permutation(len(intervals)):
                chrom, s, e = intervals[i]
                seq = extract_from_interval(chrom, s, e, jitter_rng=rng2)
                if seq is None or seq in seen:
                    continue
                seen.add(seq)
                seqs.append(seq)
                if len(seqs) >= n:
                    break
    return seqs


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    by_cat = load_ccres_by_category()
    for k, v in by_cat.items():
        print(f"{k}: {len(v):,} cCREs")

    seqs = []
    for cat in CATEGORIES:
        cat_seqs = sample_from_category(by_cat[cat], PER_CAT, rng)
        print(f"  -> {cat}: {len(cat_seqs):,} valid 200bp sequences")
        assert len(cat_seqs) == PER_CAT, f"Could not sample {PER_CAT} valid sequences for {cat}"
        seqs.extend(cat_seqs)

    rng.shuffle(seqs)  # interleave categories
    write_sequences(seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
