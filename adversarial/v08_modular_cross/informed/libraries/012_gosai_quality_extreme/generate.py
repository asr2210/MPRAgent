"""
Experiment 012 — Extreme quality filter (lfcSE<0.15).

Tests the quality asymptote. exp 005 (lfcSE<0.3) tied exp 004 (lfcSE<0.5) at
0.018 plateau — but BOTH still permit measurements with substantial noise.
lfcSE=0.15 means roughly ~30% confidence interval on a single log2FC, much
cleaner.

Method: lfcSE<0.15 + mean-activity quintile stratification.
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200
SE_THRESH = 0.15

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
        if mean_se >= SE_THRESH:
            continue
        mean_act = (k562 + hepg2 + sknsh) / 3
        records.append((seq, mean_act))

print(f"{len(records)} sequences after lfcSE<{SE_THRESH}")
if len(records) < TARGET_N:
    print(f"WARN: only {len(records)} available, need {TARGET_N}; relaxing")
    selected = records
else:
    records.sort(key=lambda x: x[1])
    n = len(records)
    per_bin = TARGET_N // 5
    selected = []
    for i in range(5):
        bin_lo = (i * n) // 5
        bin_hi = ((i + 1) * n) // 5
        bin_items = records[bin_lo:bin_hi]
        random.shuffle(bin_items)
        selected.extend(bin_items[:per_bin])

print(f"selected {len(selected)}")
seqs = [r[0] for r in selected]
random.shuffle(seqs)
assert len(seqs) == TARGET_N

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
