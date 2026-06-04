"""Exp 013 — GC-extreme DHS selection (H5 counter-test).

H5 prediction: concentration on ANY single axis hurts. GC content is the
simplest possible axis (no PWM scan needed).

Method: pre-sample 200k DHSs from non-label-aligned pool, compute GC content,
select 50k from the extreme tails (25k lowest GC + 25k highest GC). This
deliberately concentrates the library on outlier compositional DHSs.

Predicts (per H5): eval_01 < 0.6752 — both tails are non-representative.
If by some chance it lifts: H5 needs refinement (GC is a special axis).
"""
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

ROOT = Path(__file__).parent
DATA = ROOT.parent.parent / "data"
OUT = ROOT / "sequences_0.txt"

N, L = 50_000, 200
HALF = L // 2
SEED = 13
PRESCAN = 200_000

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)
print(f"  pool size: {len(df):,}")

idx = rng.choice(len(df), size=PRESCAN, replace=False)
sub = df.iloc[idx].reset_index(drop=True)

print("extracting sequences + GC...")
seqs, gc = [], []
for r in sub.itertuples(index=False):
    seq = str(fa[r.seqname][r.win_start:r.win_end]).upper()
    if "N" in seq or len(seq) != L:
        continue
    seqs.append(seq)
    gc.append((seq.count("G") + seq.count("C")) / L)
gc = np.asarray(gc)
print(f"  N-free: {len(seqs):,}  GC mean={gc.mean():.3f} std={gc.std():.3f}")
print(f"  pct GC<0.3: {(gc<0.3).mean()*100:.1f}%  pct GC>0.7: {(gc>0.7).mean()*100:.1f}%")

order = np.argsort(gc)
half = N // 2
low_idx = order[:half]    # 25k lowest GC
high_idx = order[-half:]  # 25k highest GC

selected = [seqs[i] for i in low_idx] + [seqs[i] for i in high_idx]
print(f"  low GC: {gc[low_idx].min():.3f}-{gc[low_idx].max():.3f}")
print(f"  high GC: {gc[high_idx].min():.3f}-{gc[high_idx].max():.3f}")
rng.shuffle(selected)
assert len(selected) == N

with OUT.open("w") as f:
    for s in selected:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
