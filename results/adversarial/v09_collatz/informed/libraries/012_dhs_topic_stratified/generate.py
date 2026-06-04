"""
E012: DHS topic-stratified (16 NMF components × ~3125 each).

Pivot from Gosai-centric exploration to DHS-diversity exploration.
Baseline Table 1 ranks dhs_topic/dhs_stratified above dhs_random; this
tests if my pipeline shows the same relative ordering.

Equal samples per chromatin program component forces broad regulatory
grammar coverage (cardiac, neural, etc.), which should help when eval
sets include cell types underrepresented in uniform DHS sampling.
"""
import os, gzip, numpy as np
from pyfaidx import Fasta

ROOT = "/data/users/arao/.private/MPRAgent_adversarial/runs/v09/informed_claude"
DHS = f"{ROOT}/data/DHS_Index_hg38.txt.gz"
FA = f"{ROOT}/data/hg38.fa"
SEED = 1
N_TOTAL = 50_000
LEN = 200
HALF = LEN // 2

rng = np.random.default_rng(SEED)

# Load DHS index grouped by component
print("Loading DHS index by component...", flush=True)
by_comp = {}
with gzip.open(DHS, "rt") as f:
    f.readline()
    for line in f:
        p = line.rstrip("\n").split("\t")
        comp = p[9]
        by_comp.setdefault(comp, []).append((p[0], int(p[6])))
comps = sorted(by_comp.keys())
print(f"Components: {len(comps)}", flush=True)
for c in comps:
    print(f"  {c}: {len(by_comp[c])}", flush=True)

n_per = N_TOTAL // len(comps)  # 3125 per topic
remainder = N_TOTAL - n_per * len(comps)
print(f"Sampling {n_per} per component (+{remainder} extra)", flush=True)

fa = Fasta(FA, sequence_always_upper=True, as_raw=True)
chroms = set(fa.keys())
BASES = "ACGT"

final = []
for ci, comp in enumerate(comps):
    pool = by_comp[comp]
    target = n_per + (1 if ci < remainder else 0)
    order = rng.permutation(len(pool))
    got = 0
    for ix in order:
        chrom, summit = pool[ix]
        if chrom not in chroms: continue
        s = summit - HALF; e = summit + HALF
        if s < 0 or e > len(fa[chrom]): continue
        seq = fa[chrom][s:e].upper()
        if len(seq) != LEN: continue
        nc = sum(1 for b in seq if b not in BASES)
        if nc > LEN * 0.1: continue
        if nc > 0:
            seq = "".join(b if b in BASES else BASES[rng.integers(0, 4)] for b in seq)
        final.append(seq)
        got += 1
        if got >= target: break
    print(f"  {comp}: {got}", flush=True)

# If short due to chrom/length filters, top up with random DHS
if len(final) < N_TOTAL:
    print(f"Short {N_TOTAL - len(final)}, topping up uniform", flush=True)
    all_recs = [r for v in by_comp.values() for r in v]
    order = rng.permutation(len(all_recs))
    for ix in order:
        chrom, summit = all_recs[ix]
        if chrom not in chroms: continue
        s = summit - HALF; e = summit + HALF
        if s < 0 or e > len(fa[chrom]): continue
        seq = fa[chrom][s:e].upper()
        if len(seq) != LEN: continue
        nc = sum(1 for b in seq if b not in BASES)
        if nc > LEN * 0.1: continue
        if nc > 0:
            seq = "".join(b if b in BASES else BASES[rng.integers(0, 4)] for b in seq)
        final.append(seq)
        if len(final) >= N_TOTAL: break

rng.shuffle(final)
assert len(final) == N_TOTAL, f"Got {len(final)}"

out = os.path.join(os.path.dirname(__file__), "sequences_0.txt")
with open(out, "w") as f:
    f.write("\n".join(final))
    f.write("\n")
print(f"Wrote {len(final)} to {out}", flush=True)
