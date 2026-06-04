"""
Experiment 026 — GTEX lfcSE<0.7, TOP 50K by |mean activity| (deterministic).

Removes seed variance by deterministically picking the 50K GTEX sequences with
the largest effect magnitude (within quality filter). If high-|effect|
sequences are the most informative, this should reliably beat random sampling.
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
            ["K562_log2FC", "HepG2_log2FC", "SKNSH_log2FC",
             "K562_lfcSE", "HepG2_lfcSE", "SKNSH_lfcSE"]}
    for line in f:
        parts = line.rstrip("\n").split("\t")
        if parts[proj_idx] != "GTEX":
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
        if (k_se + h_se + s_se) / 3 >= SE_CUT:
            continue
        mean_act = (k562 + hepg2 + sknsh) / 3
        records.append((seq, abs(mean_act)))

records.sort(key=lambda x: -x[1])
print(f"{len(records)} GTEX seqs after lfcSE<{SE_CUT}, top abs_mean cutoff at: {records[TARGET_N-1][1]:.3f}")
selected = records[:TARGET_N]
seqs = [r[0] for r in selected]
random.shuffle(seqs)

with open(OUT, "w") as f:
    for s in seqs:
        f.write(s + "\n")
print(f"wrote {len(seqs)} (top 50K by |mean activity|) to {OUT}")
