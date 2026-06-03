"""
E003: Top-signal DHS sites with multi-cell-type accessibility.

Pick top 50k DHS sites by (mean_signal × log(1+numsamples)). This
combines accessibility strength with cell-type breadth — sites likely
active across many cell types AND with strong open-chromatin signal.

Hypothesis: under real labels, library SNR is the dominant factor for
training quality. Top-signal DHS sites have intrinsic MPRA activity →
higher SNR → better-trained model → eval ↑ vs uniform DHS (E002=0.32).
"""
import os, gzip
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

# Load DHS index: chrom, summit, mean_signal, numsamples
print("Loading DHS index...", flush=True)
recs = []
with gzip.open(DHS_PATH, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        chrom, summit = p[0], int(p[6])
        ms = float(p[4])      # mean_signal
        ns = int(p[5])        # numsamples
        recs.append((chrom, summit, ms, ns))
print(f"Loaded {len(recs)}", flush=True)

# Compute score
score = np.array([r[2] * np.log1p(r[3]) for r in recs])
# Top 80k candidates (oversample to handle N-rich/boundary failures)
top_n = 80_000
top_idx = np.argpartition(score, -top_n)[-top_n:]
# Sort within top
top_idx = top_idx[np.argsort(-score[top_idx])]

print(f"Top {top_n} candidates by signal x log(samples)", flush=True)

fa = Fasta(FA_PATH, sequence_always_upper=True, as_raw=True)
chroms = set(fa.keys())
BASES = "ACGT"

selected = []
for ix in top_idx:
    chrom, summit, ms, ns = recs[ix]
    if chrom not in chroms:
        continue
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
    n_count = sum(1 for b in seq if b not in BASES)
    if n_count > LEN * 0.1:
        continue
    if n_count > 0:
        seq = "".join(b if b in BASES else BASES[rng.integers(0, 4)] for b in seq)
    selected.append(seq)
    if len(selected) >= N_SEQS:
        break

print(f"Selected {len(selected)}", flush=True)
assert len(selected) == N_SEQS

out_path = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out_path, "w") as f:
    f.write("\n".join(selected))
    f.write("\n")
print(f"Wrote {out_path}", flush=True)
