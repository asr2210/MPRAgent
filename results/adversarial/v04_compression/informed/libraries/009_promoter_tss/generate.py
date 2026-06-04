#!/usr/bin/env python3
"""
009_promoter_tss — 200bp windows centered on RefSeq transcript TSSes.

Motivation (Agarwal et al., Nature 2024): 200bp promoter cores
function as 'on switches' that are largely CELL-TYPE-AGNOSTIC. They
provide similar expression levels across tissues, in contrast to
enhancers which are tissue-specific. If true in this pipeline,
promoter-centered sequences should be MORE generalizable across cell
types than random genome (which mixes promoters, enhancers, repeats,
gene bodies, intergenic).

Library: 40,282 unique transcript TSSes from RefSeq on chr1–22+X+Y.
- For each TSS, draw 200bp centered on it (TSS ± 100, with random
  offset in [-50, +50] for diversity).
- Oversample: cycle through TSSes with different offsets until 50k.
- Discard windows that overlap N runs.

Generalization argument: per Agarwal 2024, promoter cores are
universal — high activity in any cell type that expresses the gene.
A library dominated by promoters should encode the universal
"transcription-machinery-engaging" grammar that any cell-type model
needs. The trade-off vs random genome: less coverage of enhancer /
intergenic / repeat content, but higher activity density per slot.

Tests:
- 009 > 002 (0.499): promoters generalize better than random genome.
  Confirms "universal regulators" hypothesis.
- 009 < 002: promoters are too narrow; random genome diversity wins.
- 009 ≈ 002: promoter cores carry no extra signal in this pipeline.

Note: this uses transcript-level TSSes (40k), oversamples with random
offsets. Many genes have multiple transcripts → multiple nearby TSSes,
which is realistic for the regulatory landscape.
"""
import gzip
from pathlib import Path

import numpy as np

SEED = 0
N = 50_000
L = 200
OFFSET = 50  # random shift ±OFFSET around TSS-centered window
DATA_DIR = Path(__file__).resolve().parents[2] / "data"
CHROMS = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]


def load_chrom(path: Path) -> str:
    with gzip.open(path, "rt") as f:
        f.readline()
        return "".join(line.strip() for line in f).upper()


def main():
    rng = np.random.default_rng(SEED)
    # Parse TSSes
    tss = []
    std = set(CHROMS)
    seen = set()
    with gzip.open(DATA_DIR / "hg38.refGene.gtf.gz", "rt") as f:
        for line in f:
            parts = line.split("\t")
            if parts[2] != "transcript":
                continue
            chrom = parts[0]
            if chrom not in std:
                continue
            start = int(parts[3])
            end = int(parts[4])
            strand = parts[6]
            t = start if strand == "+" else end
            key = (chrom, t)
            if key in seen:
                continue
            seen.add(key)
            tss.append((chrom, t))
    print(f"unique TSSes: {len(tss):,}")

    # Load all needed chromosomes
    chrom_seqs = {c: load_chrom(DATA_DIR / f"{c}.fa.gz") for c in CHROMS}
    print(f"loaded {len(chrom_seqs)} chromosomes")

    valid = set("ACGT")
    out = []
    attempts = 0
    # Shuffle TSSes for randomized iteration
    tss = list(tss)
    rng.shuffle(tss)
    i = 0
    while len(out) < N and attempts < N * 5:
        attempts += 1
        chrom, t = tss[i % len(tss)]
        i += 1
        offset = int(rng.integers(-OFFSET, OFFSET + 1))
        start = t - L // 2 + offset
        if start < 0 or start + L > len(chrom_seqs[chrom]):
            continue
        w = chrom_seqs[chrom][start : start + L]
        if not set(w).issubset(valid):
            continue
        out.append(w)
    print(f"wrote {len(out)} of {attempts} attempts "
          f"(reject {1 - len(out)/attempts:.3f})")
    assert len(out) == N
    out_path = Path(__file__).resolve().parent / "sequences_0.txt"
    with open(out_path, "w") as f:
        f.write("\n".join(out))
        f.write("\n")


if __name__ == "__main__":
    main()
