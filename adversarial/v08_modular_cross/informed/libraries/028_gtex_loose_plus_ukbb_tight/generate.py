"""
Experiment 028 — 40K GTEX-loose + 10K UKBB-tight blend.

Hypothesis: UKBB sequences need tight quality (014 vs 021 finding). A small
addition of tight-quality UKBB might add complementary signal without diluting.
40K GTEX with the loose criterion + 10K UKBB with very strict quality.
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
SEQ_LEN = 200
N_GTEX = 40_000
N_UKBB = 10_000

gtex_pool = []
ukbb_pool = []
with open(os.path.join(DATA, "Gosai_MPRA.txt")) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    proj_idx = header.index("data_project")
    cols = {c: header.index(c) for c in
            ["K562_lfcSE", "HepG2_lfcSE", "SKNSH_lfcSE"]}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        proj = parts[proj_idx]
        if proj not in ("GTEX", "UKBB"):
            continue
        seq = parts[seq_idx].upper()
        if len(seq) != SEQ_LEN or not set(seq).issubset(set("ACGT")):
            continue
        try:
            k_se = float(parts[cols["K562_lfcSE"]])
            h_se = float(parts[cols["HepG2_lfcSE"]])
            s_se = float(parts[cols["SKNSH_lfcSE"]])
        except (ValueError, IndexError):
            continue
        mean_se = (k_se + h_se + s_se) / 3
        if proj == "GTEX" and mean_se < 0.7:
            gtex_pool.append(seq)
        elif proj == "UKBB" and mean_se < 0.2:
            ukbb_pool.append(seq)

print(f"GTEX pool (lfcSE<0.7): {len(gtex_pool)}")
print(f"UKBB pool (lfcSE<0.2): {len(ukbb_pool)}")

random.shuffle(gtex_pool)
random.shuffle(ukbb_pool)
selected = gtex_pool[:N_GTEX] + ukbb_pool[:N_UKBB]
assert len(selected) == N_GTEX + N_UKBB == 50_000
random.shuffle(selected)

with open(OUT, "w") as f:
    for s in selected:
        f.write(s + "\n")
print(f"wrote {len(selected)} (40K GTEX-loose + 10K UKBB-tight) to {OUT}")
