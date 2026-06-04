"""
Experiment 013 — Signal-to-noise selection.

Select sequences with HIGH |effect|/SE in each cell type — measurements that are
both large-magnitude AND reliably measured (high t-statistic-like score).
Different from:
  - low SE alone (exp 012): selects easy/oversampled sequences
  - high |effect| alone (exp 007): selects loud sequences but mixed quality

Score: mean across cell types of |log2FC| / lfcSE.
Take top 50K by this score (no other filter), stratified across the score range
to keep diversity.
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
        # Guard against zero/tiny SE inflating score
        eps = 0.05
        snr = (abs(k562) / (k_se + eps)
             + abs(hepg2) / (h_se + eps)
             + abs(sknsh) / (s_se + eps)) / 3
        # Also include a basic noise filter to avoid garbage measurements
        mean_se = (k_se + h_se + s_se) / 3
        if mean_se >= 0.5:
            continue
        mean_act = (k562 + hepg2 + sknsh) / 3
        records.append((seq, snr, mean_act))

print(f"{len(records)} sequences after lfcSE<0.5 quality filter")

# Take top 100K by SNR, then stratify by mean_act for diversity
records.sort(key=lambda x: -x[1])
print(f"highest SNR: {records[0][1]:.2f}, snr at 50000: {records[TARGET_N-1][1]:.2f}")
print(f"snr at 100000: {records[99999][1]:.2f}")
top_pool = records[:100_000]
top_pool.sort(key=lambda x: x[2])  # by mean_act
n = len(top_pool)
per_bin = TARGET_N // 5
selected = []
for i in range(5):
    bin_lo = (i * n) // 5
    bin_hi = ((i + 1) * n) // 5
    bin_items = top_pool[bin_lo:bin_hi]
    random.shuffle(bin_items)
    selected.extend(bin_items[:per_bin])

print(f"selected {len(selected)} (top 100K by SNR, quintile-stratified by mean_act)")
seqs = [r[0] for r in selected]
random.shuffle(seqs)
assert len(seqs) == TARGET_N

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
