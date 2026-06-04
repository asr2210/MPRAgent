"""
Experiment 027 — GTEX lfcSE<0.7 with K562-activity quintile stratification.

K562 has been the dominant per-cell signal in our best libraries (consistently
2-3x the K562 score vs HepG2/SKNSH). The eval may be K562-heavy. Stratifying
by K562 activity specifically (rather than mean across cells) might match the
eval distribution better.
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200
SE_CUT = 0.7

records = []
with open(os.path.join(DATA, "Gosai_MPRA.txt")) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    proj_idx = header.index("data_project")
    cols = {c: header.index(c) for c in
            ["K562_log2FC", "K562_lfcSE", "HepG2_lfcSE", "SKNSH_lfcSE"]}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[proj_idx] != "GTEX":
            continue
        seq = parts[seq_idx].upper()
        if len(seq) != SEQ_LEN or not set(seq).issubset(set("ACGT")):
            continue
        try:
            k562 = float(parts[cols["K562_log2FC"]])
            k_se = float(parts[cols["K562_lfcSE"]])
            h_se = float(parts[cols["HepG2_lfcSE"]])
            s_se = float(parts[cols["SKNSH_lfcSE"]])
        except (ValueError, IndexError):
            continue
        if (k_se + h_se + s_se) / 3 >= SE_CUT:
            continue
        records.append((seq, k562))

print(f"{len(records)} GTEX seqs after lfcSE<{SE_CUT}")
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

print(f"selected {len(selected)} in 5 K562-activity quintiles")
seqs = [r[0] for r in selected]
random.shuffle(seqs)
assert len(seqs) == TARGET_N
with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
