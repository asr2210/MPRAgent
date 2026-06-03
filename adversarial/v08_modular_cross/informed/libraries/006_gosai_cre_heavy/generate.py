"""
Experiment 006 — CRE-heavy blend.

Composition hypothesis:
- 14K Gosai CRE (all of them, curated cCRE sub-library)
- 36K Gosai non-CRE (UKBB+GTEX), filtered lfcSE<0.3, activity-stratified per cell

If eval is CRE-enriched, this should beat exp 004/005 (~0.018).
"""
import os, random

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
DATA = os.path.join(ROOT, "data")
OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "sequences_0.txt")

random.seed(42)
TARGET_N = 50_000
SEQ_LEN = 200

cre_recs = []
other_recs = []
with open(os.path.join(DATA, "Gosai_MPRA.txt")) as f:
    header = f.readline().rstrip("\n").split("\t")
    seq_idx = header.index("sequence")
    proj_idx = header.index("data_project")
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
        proj = parts[proj_idx]
        rec = (seq, k562, hepg2, sknsh, mean_se)
        if proj == "CRE":
            # keep all CRE regardless of SE
            cre_recs.append(rec)
        else:
            if mean_se < 0.3:
                other_recs.append(rec)

print(f"CRE pool: {len(cre_recs)}")
print(f"Other pool (lfcSE<0.3): {len(other_recs)}")

# All CRE
random.shuffle(cre_recs)
selected = [r[0] for r in cre_recs]
n_cre = len(selected)
print(f"Taking all {n_cre} CRE")

# Need TARGET_N - n_cre from other pool, stratified by per-cell-type quintile
needed = TARGET_N - n_cre
print(f"Need {needed} more from other pool")

def quintiles(vals):
    s = sorted(vals)
    n = len(s)
    return [s[int(n * q / 5)] for q in (1, 2, 3, 4)]

k_b = quintiles([r[1] for r in other_recs])
h_b = quintiles([r[2] for r in other_recs])
s_b = quintiles([r[3] for r in other_recs])

def bin_idx(v, b):
    for i, bv in enumerate(b):
        if v < bv:
            return i
    return len(b)

bins = {}
for r in other_recs:
    bi = (bin_idx(r[1], k_b), bin_idx(r[2], h_b), bin_idx(r[3], s_b))
    bins.setdefault(bi, []).append(r[0])

print(f"populated {len(bins)} of 125 other-pool bins")
per_bin = needed // len(bins) + 1
other_selected = []
for bi, pool in bins.items():
    take = min(per_bin, len(pool))
    other_selected.extend(random.sample(pool, take))

print(f"From other pool: {len(other_selected)}")
random.shuffle(other_selected)
selected.extend(other_selected[:needed])

# top-up if short
all_pool = [r[0] for r in cre_recs] + [r[0] for r in other_recs]
while len(selected) < TARGET_N:
    selected.append(random.choice(all_pool))
selected = selected[:TARGET_N]
random.shuffle(selected)

with open(OUT, "w") as f:
    for s in selected:
        f.write(s + "\n")
print(f"wrote {len(selected)} to {OUT}")
