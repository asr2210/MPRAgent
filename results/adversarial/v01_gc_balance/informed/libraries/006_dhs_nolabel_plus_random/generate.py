"""Exp 006 — 25k DHS (non-label-aligned topics) + 25k uniform random.

Stacks two known-good ingredients:
  (a) exp 005: excluding label-aligned topics improved DHS performance
  (b) baseline `dhs_synth`: 50% random filler slightly improved DHS

Predicts: eval_01 ~ 0.68-0.71 if stacking works; ~0.66 if not.
"""
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

ROOT = Path(__file__).parent
DATA = ROOT.parent.parent / "data"
OUT = ROOT / "sequences_0.txt"

TOTAL, L = 50_000, 200
N_DHS, N_RND = 25_000, 25_000
HALF = L // 2
SEED = 6

EXCLUDED_TOPICS = {
    "Myeloid / erythroid", "Digestive", "Cancer / epithelial", "Neural",
}

rng = np.random.default_rng(SEED)

print("loading DHS index...")
df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                 usecols=["seqname", "summit", "component"])
df = df[~df["component"].isin(EXCLUDED_TOPICS)].reset_index(drop=True)
print(f"  pool: {len(df):,} DHSs")

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
df["win_start"] = df["summit"] - HALF
df["win_end"] = df["summit"] + HALF
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
df = df[ok].reset_index(drop=True)

# DHS half
dhs_seqs = []
attempt = 0
while len(dhs_seqs) < N_DHS and attempt < 6:
    need = N_DHS - len(dhs_seqs)
    idx = rng.choice(len(df), size=need * 2, replace=False if need * 2 <= len(df) else True)
    for i in idx:
        if len(dhs_seqs) >= N_DHS:
            break
        r = df.iloc[i]
        seq = str(fa[r.seqname][r.win_start:r.win_end]).upper()
        if "N" in seq or len(seq) != L:
            continue
        dhs_seqs.append(seq)
    attempt += 1
assert len(dhs_seqs) == N_DHS
print(f"  DHS: {len(dhs_seqs)}")

# Random half
bases = np.array(list("ACGT"))
rand_arr = bases[rng.integers(0, 4, size=(N_RND, L))]
rand_seqs = ["".join(row) for row in rand_arr]
print(f"  random: {len(rand_seqs)}")

all_seqs = dhs_seqs + rand_seqs
order = rng.permutation(len(all_seqs))
all_seqs = [all_seqs[i] for i in order]

with OUT.open('w') as f:
    for s in all_seqs:
        assert len(s) == L
        f.write(s + "\n")
print(f"wrote {len(all_seqs)} sequences x {L}bp to {OUT}")
