"""
Experiment 008 — Gosai cell-type-specific (high variance across cells).

Take 50K Gosai sequences with highest variance in activity across K562/HepG2/
SKNSH (with lfcSE<0.3 quality filter). These are the most cell-type-
discriminating sequences.

Hypothesis: if eval rewards learning cell-type-specific TF activity, cell-type-
discriminating training examples are most informative.
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
        var = statistics.pvariance([k562, hepg2, sknsh])
        records.append((seq, var))

print(f"{len(records)} sequences after quality filter")

# Sort by variance, descending, take top 50K
records.sort(key=lambda x: -x[1])
print(f"highest variance: {records[0][1]:.3f}")
print(f"50000th variance: {records[TARGET_N-1][1]:.3f}")
print(f"lowest variance: {records[-1][1]:.3f}")

seqs = [r[0] for r in records[:TARGET_N]]
random.shuffle(seqs)

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
