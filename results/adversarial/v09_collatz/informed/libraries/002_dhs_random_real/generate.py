"""
Experiment 002: DHS uniform-random baseline with real labels.

50,000 DHS sites sampled uniformly at random from the Meuleman 2020
index (~3.59M sites, hg38). Each site centered on the summit; extract
200bp window (100bp each side). Bases lowercased in soft-masked hg38
are upper-cased. If a site falls within N-rich region or near chrom
boundaries, resample. Bases outside {A,C,G,T} replaced with random.

This is a "real-labeled" counterpart to baseline `dhs_random` (which
in instructions.md is oracle-labeled at eval_01=0.7089). Tests how
much real-label noise drops DHS performance vs oracle.
"""
import os, gzip, sys
import numpy as np
from pyfaidx import Fasta

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
DHS_PATH = f"{ROOT}/data/DHS_Index_hg38.txt.gz"
FA_PATH = f"{ROOT}/data/hg38.fa"
SEED = 1
N_SEQS = 50_000
LEN = 200
HALF = LEN // 2

rng = np.random.default_rng(SEED)

# Load DHS index (chrom, summit) — skip header
print("Loading DHS index...", flush=True)
records = []
with gzip.open(DHS_PATH, "rt") as f:
    header = f.readline()
    for line in f:
        parts = line.rstrip("\n").split("\t")
        chrom, summit = parts[0], int(parts[6])
        records.append((chrom, summit))
print(f"Loaded {len(records)} DHS sites", flush=True)

# Open hg38
print("Opening hg38...", flush=True)
fa = Fasta(FA_PATH, sequence_always_upper=True, as_raw=True)
chroms = set(fa.keys())

# Filter to records with valid chrom
records = [r for r in records if r[0] in chroms]
print(f"After chrom filter: {len(records)}", flush=True)

# Sample with replacement-tolerant rejection
selected = []
attempts = 0
target = N_SEQS
indices = rng.permutation(len(records))
i = 0
BASES = "ACGT"
while len(selected) < target and i < len(indices):
    chrom, summit = records[indices[i]]
    i += 1
    attempts += 1
    start = summit - HALF
    end = summit + HALF
    if start < 0:
        continue
    try:
        chrom_len = len(fa[chrom])
    except KeyError:
        continue
    if end > chrom_len:
        continue
    seq = fa[chrom][start:end].upper()
    if len(seq) != LEN:
        continue
    # Replace non-ACGT with random
    n_count = sum(1 for b in seq if b not in BASES)
    if n_count > LEN * 0.1:  # skip if >10% N
        continue
    if n_count > 0:
        seq = "".join(b if b in BASES else BASES[rng.integers(0, 4)] for b in seq)
    selected.append(seq)

print(f"Selected {len(selected)} sequences from {attempts} attempts", flush=True)
assert len(selected) == N_SEQS, f"Got {len(selected)} not {N_SEQS}"

out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out_path, "w") as f:
    f.write("\n".join(selected))
    f.write("\n")
print(f"Wrote to {out_path}", flush=True)
