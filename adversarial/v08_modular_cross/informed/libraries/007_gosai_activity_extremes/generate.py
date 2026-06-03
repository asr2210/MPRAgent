"""
Experiment 007 — Gosai activity extremes.

Take 25K most-active + 25K most-inactive sequences (by mean activity across
cell types) from Gosai with lfcSE<0.3 filter.

Hypothesis: extremes have the clearest signal (strong motifs vs lack thereof),
which trains better motif-learning than middle-range sequences. Tests whether
information per sequence at extremes > information at typical range.
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
        mean_act = (k562 + hepg2 + sknsh) / 3
        records.append((seq, mean_act))

print(f"{len(records)} sequences after quality filter")

records.sort(key=lambda x: x[1])
n_take_low = TARGET_N // 2
n_take_high = TARGET_N - n_take_low
low = [r[0] for r in records[:n_take_low]]
high = [r[0] for r in records[-n_take_high:]]
print(f"low: {len(low)}, activity range [{records[0][1]:.2f}, {records[n_take_low-1][1]:.2f}]")
print(f"high: {len(high)}, activity range [{records[-n_take_high][1]:.2f}, {records[-1][1]:.2f}]")

selected = low + high
random.shuffle(selected)
seqs = selected[:TARGET_N]
assert len(seqs) == TARGET_N

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
