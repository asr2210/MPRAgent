#!/usr/bin/env python3
"""
002_genome_random — real human DNA, no annotation enrichment.

50k random 200bp windows uniformly sampled from a mix of hg38
chromosomes (chr1, chr17, chr19, chr22). N-containing windows are
discarded. Each chromosome contributes proportional to its callable
(non-N) length.

This isolates the contribution of "real human DNA distribution" from
"open chromatin annotation". Compare to:
- 001_synth_random (i.i.d. uniform): 0.3068 → tells us if real DNA per
  se helps.
- Later DHS-based runs: tells us if accessibility enrichment adds on
  top of real-DNA-ness.

Generalization argument: a model trained on real genomic distribution
should learn human-typical k-mer biases, repeat content, and CpG
structure. These priors transfer across cell types because they're
properties of the DNA substrate, not of any particular cellular
context. So this design should be informative regardless of which cell
types are eventually evaluated.
"""
import gzip
import os
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CHROMS = ["chr1", "chr17", "chr19", "chr22"]


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        # Skip first header line, concat the rest, uppercase
        first = f.readline()
        assert first.startswith(">"), first
        return "".join(line.strip() for line in f).upper()


def main():
    rng = np.random.default_rng(SEED)

    # Load all chromosomes
    chrom_seqs = {}
    for c in CHROMS:
        path = DATA_DIR / f"{c}.fa.gz"
        seq = load_chrom(path)
        chrom_seqs[c] = seq
        print(f"loaded {c}: {len(seq):,} bp")

    # Compute non-N lengths to weight per-chrom sampling
    total = sum(len(s) - L for s in chrom_seqs.values())
    print(f"total callable positions (length - L): {total:,}")

    out_path = Path(__file__).resolve().parent / "sequences_0.txt"

    valid = set("ACGT")
    written = 0
    attempts = 0
    target = N
    weights = np.array([len(chrom_seqs[c]) - L for c in CHROMS], dtype=np.float64)
    weights /= weights.sum()

    with open(out_path, "w") as f:
        while written < target:
            # Sample which chromosome
            n_batch = max(target - written, 1) * 2  # over-sample to absorb N rejects
            chrom_idx = rng.choice(len(CHROMS), size=n_batch, p=weights)
            for ci in chrom_idx:
                if written >= target:
                    break
                attempts += 1
                seq = chrom_seqs[CHROMS[ci]]
                start = rng.integers(0, len(seq) - L)
                window = seq[start : start + L]
                if not set(window).issubset(valid):
                    continue
                f.write(window)
                f.write("\n")
                written += 1

    print(f"wrote {written} sequences ({attempts} attempts; reject rate "
          f"{1 - written / attempts:.3f}) → {out_path}")


if __name__ == "__main__":
    main()
