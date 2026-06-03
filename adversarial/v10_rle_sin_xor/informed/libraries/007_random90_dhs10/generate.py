"""
Experiment 007: 90% random + 10% DHS-uniform.

Hypothesis: K562 saturates near 1.0 for random; DHS-only lowers K562 to ~0.92.
A small fraction of DHS may add HepG2-relevant biological content without
sacrificing K562 saturation. If the marginal HepG2 gain > marginal K562 loss,
mean_r improves slightly.
"""
import gzip, os
import numpy as np
import twobitreader

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
DHS = os.path.join(ROOT, "data", "DHS_Index_and_Vocabulary_hg38_WM20190703.txt.gz")
TWOBIT = os.path.join(ROOT, "data", "hg38.2bit")
OUT = os.path.join(HERE, "sequences_0.txt")
SEED = 1
CHROMS = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

rng = np.random.default_rng(SEED)
bases = np.array(list("ACGT"))

# 45000 random
out = []
for _ in range(45_000):
    out.append("".join(bases[rng.integers(0, 4, 200)]))

# 5000 DHS-uniform
rows = []
with gzip.open(DHS, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        if p[0] in CHROMS:
            try:
                rows.append((p[0], int(p[6])))
            except ValueError:
                pass

idxs = rng.choice(len(rows), size=int(5_000 * 1.4), replace=False)
genome = twobitreader.TwoBitFile(TWOBIT)
dhs_seqs = []
for i in idxs:
    chrom, summit = rows[i]
    start = max(0, summit - 100)
    seq = genome[chrom][start:start + 200].upper()
    if len(seq) == 200 and all(c in "ACGT" for c in seq):
        dhs_seqs.append(seq)
out.extend(dhs_seqs[:5_000])
assert len(out) == 50_000
rng.shuffle(out)

with open(OUT, "w") as f:
    for s in out:
        f.write(s + "\n")
print(f"Wrote {len(out)} sequences (45k random + 5k DHS)")
