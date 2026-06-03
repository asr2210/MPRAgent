"""Exp 030 — cCRE 8-class balanced × chrom-balanced with SMALL-chrom boost.

Exp 028 (8-class × 23-chrom balanced) = 0.6940 NEW BEST. Its biggest lift
was eval_04/09 (+0.011 vs exp 024), suggesting that under-represented
chromosomes carried distinct signal. Small/late-replicating chroms (chr13-22,
chrX) have lower gene density and different replication timing than
chr1-12 — biology that may help generalize to UNSEEN cell types.

This experiment doubles the per-cell quota for the 11 small chroms
(chr13..chr22, chrX) — they collectively get 2x weight vs chr1..chr12.

cells = 8 classes × (12 chroms × 1 + 11 chroms × 2) = 8 × 34 = 272
per_cell = ~184 → 50048 (trimmed to 50k).

If eval_04/09 keeps climbing → small-chrom signal is the right axis to push.
If it flattens/regresses → the chrom-balance lift was a one-shot, not a
continuous gradient.
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
SEED = 30

# Weight 1 for "large" chroms (chr1..chr12); weight 2 for small chroms.
LARGE = [f"chr{i}" for i in range(1, 13)]
SMALL = [f"chr{i}" for i in range(13, 23)] + ["chrX"]
CHROM_WEIGHTS = {c: 1 for c in LARGE}
CHROM_WEIGHTS.update({c: 2 for c in SMALL})
USE_CHROMS = LARGE + SMALL

rng = np.random.default_rng(SEED)

print("loading cCREs...")
df = pd.read_csv(DATA / "GRCh38-cCREs.bed", sep="\t", header=None,
                 names=["chrom", "start", "end", "rDHS_id", "cCRE_id", "class"])
df = df[df["chrom"].isin(USE_CHROMS)].reset_index(drop=True)

df["mid"] = (df["start"] + df["end"]) // 2
df["win_start"] = df["mid"] - HALF
df["win_end"] = df["mid"] + HALF

fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
chrom_len = {c: len(fa[c]) for c in df["chrom"].unique()}
ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["chrom"]], axis=1)
df = df[ok].reset_index(drop=True)

classes = sorted(df["class"].unique())
total_weight_per_class = sum(CHROM_WEIGHTS[c] for c in USE_CHROMS)
n_units = len(classes) * total_weight_per_class
per_unit = N // n_units
remainder = N - per_unit * n_units
print(f"  classes={len(classes)}, weight_units_per_class={total_weight_per_class}, "
      f"per_unit={per_unit}, remainder={remainder}")

pool = []
extras = remainder
for cls in classes:
    sub_cls = df[df["class"] == cls]
    for chrom in USE_CHROMS:
        w = CHROM_WEIGHTS[chrom]
        sub = sub_cls[sub_cls["chrom"] == chrom].reset_index(drop=True)
        quota = per_unit * w
        # distribute extras across cells
        for _ in range(w):
            if extras > 0:
                quota += 1
                extras -= 1
        if len(sub) == 0:
            sub = sub_cls.reset_index(drop=True)
        collected = []
        attempt = 0
        while len(collected) < quota and attempt < 10:
            need = quota - len(collected)
            size = min(len(sub), need * 3)
            replace = need * 3 > len(sub)
            idx = rng.choice(len(sub), size=size, replace=replace)
            for i in idx:
                if len(collected) >= quota:
                    break
                r = sub.iloc[i]
                seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
                if "N" in seq or len(seq) != L:
                    continue
                collected.append(seq)
            attempt += 1
        if len(collected) < quota:
            sub_all = sub_cls.reset_index(drop=True)
            while len(collected) < quota:
                r = sub_all.iloc[rng.integers(len(sub_all))]
                seq = str(fa[r.chrom][r.win_start:r.win_end]).upper()
                if "N" in seq or len(seq) != L:
                    continue
                collected.append(seq)
        pool.extend(collected)

assert len(pool) == N, f"got {len(pool)}"
rng.shuffle(pool)

with OUT.open("w") as f:
    for s in pool:
        f.write(s + "\n")
print(f"wrote {N} sequences x {L}bp to {OUT}")
