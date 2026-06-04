"""Exp 004 — DHS + motif-planted synthetic hybrid (50/50).

Tests H1: a multi-objective library spanning the motif-content axis beats
single-source libraries. 25k DHS (uniform) + 25k motif-planted synthetic.

Prediction: eval_07/13 will rise above pure-DHS (motif-planting helps most),
eval_04/09 should mostly retain DHS performance (DHS half still there),
eval_08 may suffer relative to a DHS+random split.
"""
import re
import numpy as np
import pandas as pd
from pathlib import Path
from pyfaidx import Fasta

ROOT = Path(__file__).parent
DATA = ROOT.parent.parent / "data"
OUT = ROOT / "sequences_0.txt"
JASPAR = DATA / "JASPAR_vert.txt"

TOTAL, L = 50_000, 200
N_DHS, N_SYN = 25_000, 25_000
HALF = L // 2
SEED = 4
MIN_MOTIFS, MAX_MOTIFS = 3, 8
rng = np.random.default_rng(SEED)


# ---------- DHS extraction ----------
def collect_dhs(n):
    df = pd.read_csv(DATA / "DHS_Index_hg38.txt.gz", sep="\t",
                     usecols=["seqname", "summit"])
    fa = Fasta(str(DATA / "hg38.fa"), as_raw=True, sequence_always_upper=True)
    chrom_len = {c: len(fa[c]) for c in df["seqname"].unique()}
    df["win_start"] = df["summit"] - HALF
    df["win_end"] = df["summit"] + HALF
    ok = (df["win_start"] >= 0) & df.apply(lambda r: r["win_end"] <= chrom_len[r["seqname"]], axis=1)
    df = df[ok].reset_index(drop=True)

    pool, attempt = [], 0
    while len(pool) < n and attempt < 6:
        need = n - len(pool)
        idx = rng.choice(len(df), size=need * 2, replace=False if need * 2 <= len(df) else True)
        for i in idx:
            if len(pool) >= n:
                break
            r = df.iloc[i]
            seq = str(fa[r.seqname][r.win_start:r.win_end]).upper()
            if "N" in seq or len(seq) != L:
                continue
            pool.append(seq)
        attempt += 1
    assert len(pool) == n
    return pool


# ---------- Motif-planted synthetic ----------
def parse_jaspar(path):
    motifs, rows, header = [], {}, None
    def flush():
        if header and len(rows) == 4:
            counts = np.stack([rows[b] for b in 'ACGT']).astype(float) + 0.01
            pwm = counts / counts.sum(axis=0, keepdims=True)
            motifs.append((header[0], header[1], pwm))
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                flush()
                parts = line[1:].split('\t')
                header = (parts[0], parts[1] if len(parts) > 1 else parts[0])
                rows = {}
            else:
                m = re.match(r'([ACGT])\s*\[(.*)\]', line)
                if m:
                    rows[m.group(1)] = np.array([float(x) for x in m.group(2).split()])
        flush()
    return motifs


COMP = str.maketrans("ACGT", "TGCA")
def revcomp(s):
    return s.translate(COMP)[::-1]


def sample_from_pwm(pwm, r):
    cdf = pwm.cumsum(axis=0)
    u = r.random(pwm.shape[1])
    idx = (u[None, :] >= cdf).sum(axis=0)
    return "".join("ACGT"[i] for i in idx)


def synth_motif_seqs(n, motifs):
    bases = np.array(list("ACGT"))
    backbone = bases[rng.integers(0, 4, size=(n, L))]
    seqs_chars = [list("".join(row)) for row in backbone]
    n_per = rng.integers(MIN_MOTIFS, MAX_MOTIFS + 1, size=n)
    motif_choices = rng.integers(0, len(motifs), size=int(n_per.sum()))
    strands = rng.random(int(n_per.sum())) < 0.5
    cur = 0
    out = []
    for i in range(n):
        seq = seqs_chars[i]
        placed = []
        for _ in range(n_per[i]):
            mi = motif_choices[cur]; strand = strands[cur]; cur += 1
            _, _, pwm = motifs[mi]
            ml = pwm.shape[1]
            if ml >= L:
                continue
            for _try in range(6):
                pos = rng.integers(0, L - ml + 1)
                end = pos + ml
                if any(not (end <= s or pos >= e) for s, e in placed):
                    continue
                inst = sample_from_pwm(pwm, rng)
                if strand:
                    inst = revcomp(inst)
                for off, ch in enumerate(inst):
                    seq[pos + off] = ch
                placed.append((pos, end))
                break
        out.append("".join(seq))
    return out


def main():
    print("collecting DHS sequences...")
    dhs_seqs = collect_dhs(N_DHS)
    print(f"  got {len(dhs_seqs)} DHS")
    print("parsing JASPAR...")
    motifs = parse_jaspar(JASPAR)
    print(f"  {len(motifs)} motifs")
    print("synthesising motif-planted sequences...")
    syn_seqs = synth_motif_seqs(N_SYN, motifs)
    print(f"  got {len(syn_seqs)} synth")
    all_seqs = dhs_seqs + syn_seqs
    # Shuffle so training sees a mix
    order = rng.permutation(len(all_seqs))
    all_seqs = [all_seqs[i] for i in order]
    with OUT.open('w') as f:
        for s in all_seqs:
            assert len(s) == L
            f.write(s + "\n")
    print(f"wrote {len(all_seqs)} sequences x {L}bp to {OUT}")


if __name__ == "__main__":
    main()
