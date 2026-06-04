#!/usr/bin/env python3
"""
011_multisource_merge — 5 sources × 10k = 50k. Maximum diversity across
regulatory element types and data modalities.

Sources:
  - cCRE PLS (promoter-like, 10k)
  - cCRE dELS (distal enhancer-like, 10k)
  - DHS summits highest mean_signal (10k from top of Meuleman)
  - Sharpr-MPRA top 10k by mean activity (real measured strong enhancers)
  - Sharpr-MPRA bottom 10k by mean activity (real measured silencers)

Rationale: covers four orthogonal axes — promoter vs enhancer (cCRE classes),
accessibility-prioritized (DHS), and *measured* activity poles (Sharpr).
A model trained on this should be exposed to virtually every regulatory
sequence flavor relevant to a generalizing MPRA predictor.
"""
import gzip
import random
from pathlib import Path

import h5py
import numpy as np
from pyfaidx import Fasta

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
DATA = ROOT / "data"
CCRE = DATA / "GRCh38-cCREs.bed"
DHS = DATA / "DHS_Index_hg38.txt.gz"
SHARPR = DATA / "sharpr_train.hdf5"
GENOME = DATA / "hg38.fa"
OUT = HERE / "sequences_0.txt"

L = 200
PER_SRC = 10_000
SEED = 0
ACGT = set("ACGT")
ACGT_STR = "ACGT"


def base_class(label):
    for p in label.split(","):
        if p in ("PLS", "pELS", "dELS", "CTCF-only", "DNase-H3K4me3"):
            return p
    return None


def extract_genomic(coords, fa, n, label):
    out = []
    for chrom, mid in coords:
        if len(out) >= n:
            break
        s = mid - L // 2
        e = s + L
        if s < 0:
            continue
        try:
            clen = len(fa[chrom])
        except KeyError:
            continue
        if e > clen:
            continue
        seq = str(fa[chrom][s:e]).upper()
        if len(seq) != L or not set(seq).issubset(ACGT):
            continue
        out.append(seq)
    print(f"  {label}: {len(out)}/{n}")
    return out


def load_ccre(target_cls, rng):
    candidates = []
    with open(CCRE) as f:
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            if "_" in chrom or chrom == "chrM":
                continue
            if base_class(parts[-1]) != target_cls:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            candidates.append((chrom, (start + end) // 2))
    rng.shuffle(candidates)
    return candidates


def load_dhs_top(n_pool, rng):
    rows = []
    with gzip.open(DHS, "rt") as f:
        next(f)
        for line in f:
            parts = line.rstrip().split("\t")
            chrom = parts[0]
            if "_" in chrom or chrom == "chrM":
                continue
            try:
                summit = int(parts[6])
                signal = float(parts[4])
            except ValueError:
                continue
            rows.append((signal, chrom, summit))
    rows.sort(key=lambda r: r[0], reverse=True)
    return [(c, m) for _, c, m in rows[:n_pool]]


def load_sharpr(rng):
    with h5py.File(SHARPR, "r") as f:
        Y = f["Y/output"][:]
        mean_act = Y.mean(axis=1)
        order = np.argsort(mean_act)
        bot_idx = order[:PER_SRC * 2]
        top_idx = order[-PER_SRC * 2:]
        all_idx = np.concatenate([bot_idx, top_idx])
        all_idx_sorted = np.sort(all_idx)
        X_sel = f["X/sequence"][all_idx_sorted]
    # remap original to sorted positions
    pos = {orig: i for i, orig in enumerate(all_idx_sorted)}
    top_in_sel = [pos[i] for i in top_idx]
    bot_in_sel = [pos[i] for i in bot_idx]

    idx = X_sel.argmax(axis=-1)
    has_base = X_sel.sum(axis=-1) > 0

    def decode_pad(positions, n, label):
        out = []
        for p in positions:
            if len(out) >= n:
                break
            row = idx[p]
            mask = has_base[p]
            chars = [ACGT_STR[row[j]] if mask[j] else rng.choice(ACGT_STR) for j in range(145)]
            core = "".join(chars)
            left = "".join(rng.choices(ACGT_STR, k=27))
            right = "".join(rng.choices(ACGT_STR, k=28))
            out.append(left + core + right)
        print(f"  Sharpr-{label}: {len(out)}/{n}")
        return out

    return decode_pad(top_in_sel, PER_SRC, "top"), decode_pad(bot_in_sel, PER_SRC, "bot")


def main():
    rng = random.Random(SEED)
    fa = Fasta(str(GENOME), as_raw=True, sequence_always_upper=True)

    print("Loading cCRE-PLS...")
    pls = load_ccre("PLS", rng)
    print("Loading cCRE-dELS...")
    dels_ = load_ccre("dELS", rng)
    print("Loading DHS top...")
    dhs = load_dhs_top(PER_SRC * 3, rng)

    out = []
    out.extend(extract_genomic(pls, fa, PER_SRC, "cCRE-PLS"))
    out.extend(extract_genomic(dels_, fa, PER_SRC, "cCRE-dELS"))
    out.extend(extract_genomic(dhs, fa, PER_SRC, "DHS-top"))

    print("Loading Sharpr...")
    s_top, s_bot = load_sharpr(rng)
    out.extend(s_top)
    out.extend(s_bot)

    assert len(out) == PER_SRC * 5, f"got {len(out)}"
    rng.shuffle(out)
    with open(OUT, "w") as f:
        for s in out:
            f.write(s + "\n")
    print(f"Wrote {len(out)} -> {OUT}")


if __name__ == "__main__":
    main()
