"""
Experiment 005 — Gosai with tighter quality + per-cell-type stratification.

Push the quality+diversity dial:
- Restrict to mean lfcSE < 0.3 (stricter than 0.5 in exp 004)
- Stratify by per-cell-type activity quintile: 5x5x5 = 125 strata in 3D,
  but we use 5 bins per cell type independently (15 strata total).
- Sample uniformly within each cell-type quintile.

Hypothesis: stricter quality + cell-type-specific stratification gives cleaner
labels and better coverage of cell-type-specific activity space.
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200

records = []
with open(os.path.join(DATA, "Gosai_MPRA.txt")) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    cols = {c: header.index(c) for c in
            ["K562_log2FC", "HepG2_log2FC", "SKNSH_log2FC",
             "K562_lfcSE", "HepG2_lfcSE", "SKNSH_lfcSE"]}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        seq = parts[seq_idx].upper()
        if len(seq) != SEQ_LEN or not set(seq).issubset(set("ACGT")):
            continue
        try:
            k562 = float(parts[cols["K562_log2FC"]])
            hepg2 = float(parts[cols["HepG2_log2FC"]])
            sknsh = float(parts[cols["SKNSH_log2FC"]])
            k_se = float(parts[cols["K562_lfcSE"]])
            h_se = float(parts[cols["HepG2_lfcSE"]])
            s_se = float(parts[cols["SKNSH_lfcSE"]])
        except (ValueError, IndexError):
            continue
        mean_se = (k_se + h_se + s_se) / 3
        if mean_se >= 0.3:
            continue
        records.append((seq, k562, hepg2, sknsh))

print(f"{len(records)} sequences after lfcSE < 0.3 filter")

# Determine quintile boundaries per cell type
def quintiles(vals):
    s = sorted(vals)
    n = len(s)
    return [s[int(n * q / 5)] for q in (1, 2, 3, 4)]

k_boundaries = quintiles([r[1] for r in records])
h_boundaries = quintiles([r[2] for r in records])
s_boundaries = quintiles([r[3] for r in records])
print(f"K562 quintile boundaries: {k_boundaries}")
print(f"HepG2 quintile boundaries: {h_boundaries}")
print(f"SKNSH quintile boundaries: {s_boundaries}")

def bin_idx(v, boundaries):
    for i, b in enumerate(boundaries):
        if v < b:
            return i
    return len(boundaries)

# Bin each record into 3D (5,5,5) = 125 cells
bins = {}
for rec in records:
    seq, k, h, s = rec
    bi = (bin_idx(k, k_boundaries), bin_idx(h, h_boundaries), bin_idx(s, s_boundaries))
    bins.setdefault(bi, []).append(seq)

print(f"populated {len(bins)} of 125 bins")
# Aim for ~400 per bin to reach 50000, but some bins will be smaller (off-diagonal)
target_per_bin = TARGET_N // 125 + 1  # ~401
selected = []
for bi, pool in bins.items():
    take = min(target_per_bin, len(pool))
    selected.extend(random.sample(pool, take))

print(f"after stratified sampling: {len(selected)}")

# Top-up randomly if short
all_seqs_pool = [r[0] for r in records]
random.shuffle(all_seqs_pool)
pi = 0
while len(selected) < TARGET_N:
    selected.append(all_seqs_pool[pi])
    pi += 1

# Trim if too many
random.shuffle(selected)
seqs = selected[:TARGET_N]

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
