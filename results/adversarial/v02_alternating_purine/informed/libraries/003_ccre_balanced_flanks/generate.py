"""003_ccre_balanced_flanks — cCRE 5-class balanced + paired genomic flanks.

Hypothesis: A library that gives the model EQUAL exposure to each of the 5
canonical ENCODE cCRE classes (PLS, pELS, dELS, CTCF-only, DNase-H3K4me3) +
paired flanks should reach blind's best (~0.17). Clean class-balance vs.
blind 013's "asymmetric" hack (15k uniform + 5k CTCF + 5k DNH3).

Design:
- 5,000 positives per cCRE main class * 5 classes = 25,000 positives
- 25,000 paired flanks (±1500-3000bp from each positive center, no cCRE overlap)
- Total 50,000

Rationale:
- Class balance gives the model equal information per regulatory grammar.
- Promoter (PLS) is starved in dELS-dominated natural distribution. Equal
  share should help promoter-heavy evals (blind found 013 lost eval_13
  because it starved PLS).
- Paired flanks supply informative negatives (Avg ~50% GC contrast).
"""
from __future__ import annotations

import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
sys.path.insert(0, ROOT)
from data_utils import _fasta  # noqa: E402

DATA = os.path.join(ROOT, 'data')
N_PER_CLASS = 5_000
N_POS = N_PER_CLASS * 5
N_FLANK = 25_000
WINDOW = 200
HALF = WINDOW // 2
SEED = 0
FLANK_MIN = 1500
FLANK_MAX = 3000
CLASSES = ['PLS', 'pELS', 'dELS', 'CTCF-only', 'DNase-H3K4me3']


def load_ccres() -> pd.DataFrame:
    df = pd.read_csv(os.path.join(DATA, 'GRCh38-cCREs.bed'), sep='\t',
                     header=None,
                     names=['chrom', 'start', 'end', 'rDHS', 'accession', 'classes'])
    df['main_class'] = df['classes'].str.split(',').str[0]
    df['center'] = ((df['start'] + df['end']) // 2).astype(int)
    return df


def build_overlap_lookup(df: pd.DataFrame):
    """All cCRE intervals indexed per-chrom for flank rejection."""
    out = {}
    for chrom, g in df.groupby('chrom', sort=False):
        starts = g['start'].to_numpy()
        ends = g['end'].to_numpy()
        order = np.argsort(starts)
        out[chrom] = (starts[order], ends[order])
    return out


def overlaps(chrom, ws, we, lookup):
    if chrom not in lookup:
        return False
    starts, ends = lookup[chrom]
    j = np.searchsorted(starts, we, side='left')
    i = j - 1
    while i >= 0 and starts[i] > ws - 5000:
        if ends[i] > ws and starts[i] < we:
            return True
        i -= 1
    return False


def sample_positives(df: pd.DataFrame, n_per_class: int, seed: int):
    """Sample n_per_class centers per class. Returns DataFrame and seqs."""
    rng = np.random.default_rng(seed)
    fa = _fasta()
    out_rows, out_seqs = [], []
    n_n = n_short = n_dup = 0
    seen = set()
    for klass in CLASSES:
        pool = df[df['main_class'] == klass]
        n_need = n_per_class
        oversample = int(n_need * 1.3)
        idx = rng.choice(len(pool), size=min(oversample, len(pool)), replace=False)
        kept = 0
        for i in idx:
            if kept >= n_need:
                break
            row = pool.iloc[int(i)]
            ws, we = row['center'] - HALF, row['center'] + HALF
            try:
                s = str(fa[row['chrom']][ws:we]).upper()
            except Exception:
                n_short += 1
                continue
            if len(s) != WINDOW:
                n_short += 1
                continue
            if 'N' in s:
                n_n += 1
                continue
            if s in seen:
                n_dup += 1
                continue
            seen.add(s)
            out_seqs.append(s)
            out_rows.append((row['chrom'], int(row['center'])))
            kept += 1
        print(f'  {klass}: {kept:,}', flush=True)
    print(f'  totals: pos={len(out_seqs):,} N={n_n} short={n_short} dup={n_dup}',
          flush=True)
    return out_rows, out_seqs


def sample_flanks(positives, lookup, fa, seed, n_target):
    rng = np.random.default_rng(seed + 7)
    out, seen = [], set()
    n_overlap = n_n = n_dup = n_short = 0
    idx = rng.permutation(len(positives))
    pass_no = 0
    t0 = time.time()
    while len(out) < n_target and pass_no < 6:
        pass_no += 1
        for i in idx:
            if len(out) >= n_target:
                break
            chrom, center = positives[int(i)]
            for _ in range(4):
                sign = 1 if rng.random() < 0.5 else -1
                off = sign * int(rng.integers(FLANK_MIN, FLANK_MAX + 1))
                c2 = center + off
                ws, we = c2 - HALF, c2 + HALF
                if ws < 0:
                    n_short += 1
                    continue
                if overlaps(chrom, ws, we, lookup):
                    n_overlap += 1
                    continue
                try:
                    s = str(fa[chrom][ws:we]).upper()
                except Exception:
                    n_short += 1
                    continue
                if len(s) != WINDOW:
                    n_short += 1
                    continue
                if 'N' in s:
                    n_n += 1
                    continue
                if s in seen:
                    n_dup += 1
                    continue
                seen.add(s)
                out.append(s)
                break
        idx = rng.permutation(len(positives))
    print(f'  flanks kept {len(out):,}/{n_target} overlap={n_overlap} '
          f'N={n_n} dup={n_dup} short={n_short} in {time.time()-t0:.1f}s',
          flush=True)
    if len(out) < n_target:
        raise RuntimeError(f'only {len(out)} flanks accepted')
    return out


def main():
    t0 = time.time()
    print('loading cCREs...', flush=True)
    ccres = load_ccres()
    print(f'  {len(ccres):,} cCREs', flush=True)

    print(f'sampling {N_POS:,} positives ({N_PER_CLASS:,}/class)...', flush=True)
    pos_meta, pos_seqs = sample_positives(ccres, N_PER_CLASS, SEED)
    assert len(pos_seqs) == N_POS

    print('building cCRE overlap lookup...', flush=True)
    lookup = build_overlap_lookup(ccres)

    print(f'sampling {N_FLANK:,} paired flanks...', flush=True)
    fa = _fasta()
    flank_seqs = sample_flanks(pos_meta, lookup, fa, SEED, N_FLANK)

    all_seqs = pos_seqs + flank_seqs
    assert len(all_seqs) == 50_000
    rng = np.random.default_rng(SEED + 99)
    rng.shuffle(all_seqs)

    out_path = os.path.join(HERE, 'sequences_0.txt')
    with open(out_path, 'w') as f:
        f.write('\n'.join(all_seqs) + '\n')
    print(f'wrote {out_path} in {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
