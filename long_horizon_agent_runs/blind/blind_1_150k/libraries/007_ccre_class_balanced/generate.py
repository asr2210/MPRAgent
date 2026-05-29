"""Experiment 007 — Class-balanced cCREs.

Forces equal representation of four cCRE class buckets:
  1. PLS-types     — PLS, PLS+CTCF
  2. pELS-types    — pELS, pELS+CTCF
  3. dELS-types    — dELS, dELS+CTCF
  4. Other         — CTCF-only, DNase-H3K4me3 with or without CTCF

Natural cCRE distribution is ~74% dELS; balancing exposes the model to
typically underrepresented promoter (PLS) and insulator (CTCF-only)
grammars. Hypothesis: this lifts hard evals 11/12 if they test promoter
or insulator function.

37,500 per bucket = 150,000 total. If a bucket has fewer survivors after
filtering, sample with replacement to fill.
"""
from __future__ import annotations

import os
import sys

import numpy as np

N_SEQS = 150_000
PER_BUCKET = N_SEQS // 4  # 37,500
SEQ_LEN = 200
HALF = SEQ_LEN // 2
SEED = 6
DATA_DIR = "/data/users/arao/mpra_autoresearch/data"
BED_PATH = os.path.join(DATA_DIR, "encode_ccres_hg38.bed")
OUT_PATH = os.path.join(os.path.dirname(__file__), "sequences.txt")

PLS_CLASSES = {"PLS", "PLS,CTCF-bound"}
PELS_CLASSES = {"pELS", "pELS,CTCF-bound"}
DELS_CLASSES = {"dELS", "dELS,CTCF-bound"}
OTHER_CLASSES = {"CTCF-only,CTCF-bound", "DNase-H3K4me3", "DNase-H3K4me3,CTCF-bound"}

BUCKETS = [
    ("PLS-types", PLS_CLASSES),
    ("pELS-types", PELS_CLASSES),
    ("dELS-types", DELS_CLASSES),
    ("Other", OTHER_CLASSES),
]


def read_fasta_seq(path: str) -> str:
    chunks: list[str] = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                if chunks:
                    break
                continue
            chunks.append(line.rstrip())
    return "".join(chunks)


def main() -> None:
    # Load chromosome sequences (lazy: build set of chroms first).
    bucket_cres: dict[str, list[tuple[str, int]]] = {name: [] for name, _ in BUCKETS}
    chroms_needed: set[str] = set()
    bucket_for_class = {}
    for name, classes in BUCKETS:
        for c in classes:
            bucket_for_class[c] = name

    with open(BED_PATH) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            mid = (int(parts[1]) + int(parts[2])) // 2
            cls = parts[5]
            bucket = bucket_for_class.get(cls)
            if bucket is None:
                continue  # ignore unknown / cross-class
            bucket_cres[bucket].append((chrom, mid))
            chroms_needed.add(chrom)

    for name, _ in BUCKETS:
        print(f"  bucket {name}: {len(bucket_cres[name]):,} cCREs", file=sys.stderr)

    chrom_seq: dict[str, str] = {}
    for c in sorted(chroms_needed):
        path = os.path.join(DATA_DIR, f"{c}.fa")
        if os.path.exists(path):
            chrom_seq[c] = read_fasta_seq(path)

    valid = set("ACGTacgt")

    def keep(chrom: str, mid: int) -> str | None:
        seq = chrom_seq.get(chrom)
        if seq is None:
            return None
        lo, hi = mid - HALF, mid + HALF
        if lo < 0 or hi > len(seq):
            return None
        w = seq[lo:hi]
        if any(c not in valid for c in w):
            return None
        if sum(1 for c in w if c.islower()) > SEQ_LEN // 2:
            return None
        return w.upper()

    rng = np.random.default_rng(SEED)
    all_sample: list[str] = []
    for name, _ in BUCKETS:
        pool = [s for chrom, mid in bucket_cres[name] if (s := keep(chrom, mid)) is not None]
        print(f"  bucket {name} survivors: {len(pool):,}", file=sys.stderr)
        if len(pool) >= PER_BUCKET:
            idx = rng.choice(len(pool), size=PER_BUCKET, replace=False)
        else:
            print(
                f"    WARN: bucket {name} short ({len(pool)} < {PER_BUCKET}); "
                "sampling with replacement",
                file=sys.stderr,
            )
            idx = rng.choice(len(pool), size=PER_BUCKET, replace=True)
        sample = [pool[i] for i in idx]
        all_sample.extend(sample)

    assert len(all_sample) == N_SEQS
    rng.shuffle(all_sample)
    assert all(len(s) == SEQ_LEN and set(s) <= set("ACGT") for s in all_sample[:10])

    with open(OUT_PATH, "w") as f:
        f.write("\n".join(all_sample))
        f.write("\n")

    print(f"Wrote {N_SEQS:,} sequences to {OUT_PATH}", file=sys.stderr)


if __name__ == "__main__":
    main()
