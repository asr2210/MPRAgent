"""
E009: 50/50 mixture — 25K Gosai random + 25K DHS uniform.

Tests if combining two complementary data sources improves on either alone.
Gosai = pre-tested MPRA sequences (proven activity). DHS = chromatin
accessibility sites (broad regulatory grammar). Either alone gives ~0.32;
hypothesis is that combining them gives small synergy from regulatory
diversity that neither source alone covers.

Single seed for both halves (deterministic).
"""
import os, gzip, numpy as np
from pyfaidx import Fasta

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
GOSAI = f"{ROOT}/data/gosai_mpra.txt"
DHS = f"{ROOT}/data/DHS_Index_hg38.txt.gz"
FA = f"{ROOT}/data/hg38.fa"
SEED = 1
N_HALF = 25_000
LEN = 200
HALF = LEN // 2

rng = np.random.default_rng(SEED)

# Gosai half
print("Loading Gosai...", flush=True)
g = []
with open(GOSAI) as f:
    f.readline()
    for line in f:
        c = line.rstrip("\n").split("\t")
        if len(c) < 12 or len(c[11]) != 200: continue
        if any(b not in "ACGT" for b in c[11]): continue
        g.append(c[11])
print(f"Gosai pool: {len(g)}", flush=True)
g_sel = [g[i] for i in rng.permutation(len(g))[:N_HALF]]

# DHS half
print("Loading DHS index...", flush=True)
recs = []
with gzip.open(DHS, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        recs.append((p[0], int(p[6])))
print(f"DHS sites: {len(recs)}", flush=True)
fa = Fasta(FA, sequence_always_upper=True, as_raw=True)
chroms = set(fa.keys())
BASES = "ACGT"
order = rng.permutation(len(recs))
d_sel = []
for ix in order:
    chrom, summit = recs[ix]
    if chrom not in chroms: continue
    s = summit - HALF; e = summit + HALF
    if s < 0 or e > len(fa[chrom]): continue
    seq = fa[chrom][s:e].upper()
    if len(seq) != LEN: continue
    nc = sum(1 for b in seq if b not in BASES)
    if nc > LEN*0.1: continue
    if nc > 0:
        seq = "".join(b if b in BASES else BASES[rng.integers(0,4)] for b in seq)
    d_sel.append(seq)
    if len(d_sel) >= N_HALF: break

print(f"Gosai={len(g_sel)} DHS={len(d_sel)}", flush=True)

final = g_sel + d_sel
rng.shuffle(final)
assert len(final) == 50_000

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {out}", flush=True)
