"""
Experiment 009 — Mixed: extremes + stratified middle.

Combines:
  - 25K activity extremes (12.5K top + 12.5K bottom by mean activity, lfcSE<0.3)
  - 25K stratified middle (5 quintiles within the middle 60% of activity)

Hypothesis: extremes boost eval_04/09 (per exp 007), stratified middle preserves
eval_01 (per exp 005). A blended library should improve mean correlation across
all 14 eval sets.
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200
N_EXTREMES = 25_000  # 12.5K top + 12.5K bottom
N_MIDDLE = 25_000    # stratified within the middle range

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
        mean_act = (k562 + hepg2 + sknsh) / 3
        records.append((seq, mean_act))

print(f"{len(records)} sequences after quality filter")
records.sort(key=lambda x: x[1])
n = len(records)

# Pick bottom 12.5K and top 12.5K as extremes
n_bot = N_EXTREMES // 2
n_top = N_EXTREMES - n_bot
extremes_bot = records[:n_bot]
extremes_top = records[-n_top:]
extremes = extremes_bot + extremes_top

# Middle 60%: indices [0.2n, 0.8n). Stratify into 5 quintiles within this band.
mid_start = int(0.2 * n)
mid_end = int(0.8 * n)
middle_pool = records[mid_start:mid_end]
random.shuffle(middle_pool)
# Resort middle by activity for stratification
middle_pool.sort(key=lambda x: x[1])
m = len(middle_pool)
per_bin = N_MIDDLE // 5
middle_seqs = []
for i in range(5):
    bin_lo = (i * m) // 5
    bin_hi = ((i + 1) * m) // 5
    bin_items = middle_pool[bin_lo:bin_hi]
    random.shuffle(bin_items)
    middle_seqs.extend(bin_items[:per_bin])

print(f"extremes: {len(extremes)} (low activity range "
      f"[{extremes_bot[0][1]:.2f}, {extremes_bot[-1][1]:.2f}], high "
      f"[{extremes_top[0][1]:.2f}, {extremes_top[-1][1]:.2f}])")
print(f"middle: {len(middle_seqs)} (activity range "
      f"[{middle_pool[0][1]:.2f}, {middle_pool[-1][1]:.2f}])")

selected = [r[0] for r in extremes + middle_seqs]
assert len(selected) == TARGET_N, f"got {len(selected)} want {TARGET_N}"
random.shuffle(selected)

with open(OUT, "w") as f:
    for s in selected:
        f.write(s + "\n")
print(f"wrote {len(selected)} to {OUT}")
