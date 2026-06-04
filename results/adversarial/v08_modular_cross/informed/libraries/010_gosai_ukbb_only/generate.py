"""
Experiment 010 — Gosai UKBB sub-source only, quality + stratified.

Gosai composition: 14K CRE + 446K GTEX + 338K UKBB. Random Gosai sampling picks
~58% GTEX + 42% UKBB. This experiment isolates UKBB (GWAS variant-centric MPRA)
to see if sub-source matters for eval.

Method:
  - data_project == "UKBB" only
  - lfcSE<0.3 quality filter
  - Stratify by mean activity quintile, 10K per bin
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
    proj_idx = header.index("data_project")
    cols = {c: header.index(c) for c in
            ["K562_log2FC", "HepG2_log2FC", "SKNSH_log2FC",
             "K562_lfcSE", "HepG2_lfcSE", "SKNSH_lfcSE"]}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[proj_idx] != "UKBB":
            continue
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

print(f"{len(records)} UKBB sequences after lfcSE<0.3")
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

print(f"selected {len(selected)} in 5 quintile bins")
seqs = [r[0] for r in selected]
random.shuffle(seqs)
assert len(seqs) == TARGET_N

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} to {OUT}")
