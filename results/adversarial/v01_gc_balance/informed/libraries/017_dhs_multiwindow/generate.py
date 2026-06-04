"""Exp 017 — multi-window per DHS source (paired views).

H7: uniform within filter wins. Augmentation that breaks an invariance the
model lacks would help. Exp 009 showed positional jitter (random offset per
sequence) is flat — but that test gives ONE view per DHS at a random offset,
not paired views of the same source.

This experiment: sample 25k unique DHSs from non-label-aligned pool, emit
TWO windows per DHS at offsets (-50, +50). Each DHS contributes 2 sequences
that overlap by 100bp around the summit but otherwise differ. Library
= 50k from 25k biological sources, paired views.

If eval_01 > 0.6752 → paired views per source help (model learns to invariance
better when it sees the same regulatory content shifted).
If flat → biology is the bottleneck, not augmentation density.
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
OFFSET = 50
N_DHS = N // 2  # 25k DHSs
SEED = 17

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
df["min_start"] = df["summit"] - HALF - OFFSET
df["max_end"] = df["summit"] + HALF + OFFSET
ok = (df["min_start"] >= 0) & df.apply(lambda r: r["max_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)
print(f"  pool after edge filter: {len(df):,}")

pool, attempt = [], 0
while len(pool) < N and attempt < 6:
    need_dhs = (N - len(pool) + 1) // 2
    idx = rng.choice(len(df), size=need_dhs * 2, replace=False if need_dhs * 2 <= len(df) else True)
    for i in idx:
        if len(pool) >= N:
            break
        r = df.iloc[i]
        # left window (offset -OFFSET) and right window (offset +OFFSET)
        for off in (-OFFSET, +OFFSET):
            center = int(r.summit) + off
            s, e = center - HALF, center + HALF
            seq = str(fa[r.seqname][s:e]).upper()
            if "N" in seq or len(seq) != L:
                continue
            pool.append(seq)
            if len(pool) >= N:
                break
    attempt += 1
assert len(pool) == N

rng.shuffle(pool)  # interleave the paired views

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
