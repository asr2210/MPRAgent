"""
Experiment 004 — Gosai with quality filter + activity stratification.

Strategy:
1. Use ALL Gosai sequences (any chromosome) — the chr-holdout filter showed no
   benefit in exp 003, so it likely doesn't matter.
2. Filter to 200bp ACGT and mean lfcSE < 0.5 (high-confidence measurements).
3. Stratify by activity quintiles within each cell type to ensure full dynamic
   range coverage. Use mean activity across cell types for stratification simplicity.

If eval_01 jumps to >0.05, quality filtering helps. If it stays near 0, the
problem isn't sampling quality.
"""
import os, random
import statistics

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
    k562_idx = header.index("K562_log2FC")
    hepg2_idx = header.index("HepG2_log2FC")
    sknsh_idx = header.index("SKNSH_log2FC")
    k562_se = header.index("K562_lfcSE")
    hepg2_se = header.index("HepG2_lfcSE")
    sknsh_se = header.index("SKNSH_lfcSE")
    for line in f:
        parts = line.rstrip("\n").split("\t")
        seq = parts[seq_idx].upper()
        if len(seq) != SEQ_LEN or not set(seq).issubset(set("ACGT")):
            continue
        try:
            k562 = float(parts[k562_idx])
            hepg2 = float(parts[hepg2_idx])
            sknsh = float(parts[sknsh_idx])
            k_se = float(parts[k562_se])
            h_se = float(parts[hepg2_se])
            s_se = float(parts[sknsh_se])
        except (ValueError, IndexError):
            continue
        mean_se = (k_se + h_se + s_se) / 3
        if mean_se > 0.5:
            continue
        mean_act = (k562 + hepg2 + sknsh) / 3
        records.append((seq, mean_act))

print(f"{len(records)} sequences after quality filter")

# Stratify by activity quintiles
records.sort(key=lambda x: x[1])
n = len(records)
bin_size = n // 5
target_per_bin = TARGET_N // 5

selected = []
for b in range(5):
    start = b * bin_size
    end = (b + 1) * bin_size if b < 4 else n
    bin_pool = records[start:end]
    chosen = random.sample(bin_pool, min(target_per_bin, len(bin_pool)))
    selected.extend(chosen)
    print(f"bin {b}: pool {len(bin_pool)}, took {len(chosen)}, "
          f"activity range [{records[start][1]:.2f}, {records[end-1][1]:.2f}]")

# Fill remainder if any
while len(selected) < TARGET_N:
    selected.append(random.choice(records))

random.shuffle(selected)
seqs = [s for s, _ in selected[:TARGET_N]]

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
