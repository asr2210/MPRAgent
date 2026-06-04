"""007_multisource_diverse — eval-targeted multi-source library.

Hypothesis (Theory v4): No single source covers all 14 evals. Combining
sources targeted at each eval should raise mean_r above any single-source
recipe. Per eval_mapping.json:
- eval_01/04: chr-held-out genomic ground truth -> cCRE + DHS coverage
- eval_02/03/05/06/11/12: UKBB/GTEx variant MPRA -> Table_S2 + cCRE
- eval_07: SEI chr7,13 -> Table_S2 (helped previously)
- eval_08: synthetic -> random ACGT
- eval_10: dhs_chr7_13 -> DHS-topic
- eval_13: genomic_chr7_13 -> mixed genomic
- eval_14/09: oracle variants of above

Design (all from non-eval chroms where possible):
- 10k cCRE uniform (multi-class enhancer/promoter coverage)
- 4k cCRE CTCF-only boost (insulator grammar, K562 differentiation)
- 4k cCRE DNase-H3K4me3 boost (proximal regulatory)
- 5k Table_S2 (UKBB+GTEx variant MPRA)
- 4k DHS-topic (cell-type-specific accessibility)
- 3k random synthetic ACGT (targets eval_08)
- 20k paired flanks (one per genomically-locatable positive)
Total = 50k
"""
from __future__ import annotations

import gzip
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
WINDOW = 200
HALF = WINDOW // 2
SEED = 0
FLANK_MIN = 1500
FLANK_MAX = 3000

# Eval chromosomes
EVAL_CHROMS = {'chr7', 'chr13', 'chr19', 'chr21', 'chrX'}
EVAL_CHROMS_BARE = {'7', '13', '19', '21', 'X'}
MPRA_PATH = '/data/users/arao/.private/mpra_exp/data/Table_S2__MPRA_dataset.txt'

N_CCRE_UNIFORM = 10_000
N_CCRE_CTCF = 4_000
N_CCRE_DNH3 = 4_000
N_TABLE_S2 = 5_000
N_DHS_TOPIC = 4_000
N_SYNTH = 3_000
N_FLANK = 20_000


def load_ccres():
    df = pd.read_csv(os.path.join(DATA, 'GRCh38-cCREs.bed'), sep='\t',
                     header=None,
                     names=['chrom', 'start', 'end', 'rDHS', 'accession', 'classes'])
    df['main_class'] = df['classes'].str.split(',').str[0]
    df['center'] = ((df['start'] + df['end']) // 2).astype(int)
    return df


def build_overlap_lookup(df, chrom_col):
    out = {}
    for chrom, g in df.groupby(chrom_col, sort=False):
        starts = g['start'].to_numpy() if 'start' in g.columns else None
        ends = g['end'].to_numpy() if 'end' in g.columns else None
        if starts is None:
            continue
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


def extract_centers(df, indices, fa, seen):
    out_seqs, out_meta = [], []
    for i in indices:
        row = df.iloc[int(i)]
        ws, we = row['center'] - HALF, row['center'] + HALF
        try:
            s = str(fa[row['chrom']][ws:we]).upper()
        except Exception:
            continue
        if len(s) != WINDOW or 'N' in s or s in seen:
            continue
        seen.add(s)
        out_seqs.append(s)
        out_meta.append((row['chrom'], int(row['center'])))
    return out_meta, out_seqs


def sample_cCRE(df, seed, seen):
    rng = np.random.default_rng(seed)
    fa = _fasta()
    print('  cCRE uniform 10k...', flush=True)
    pool_u = df.reset_index(drop=True)
    idx = rng.choice(len(pool_u), size=int(N_CCRE_UNIFORM * 1.3), replace=False)
    mu, su = extract_centers(pool_u, idx, fa, seen)
    mu = mu[:N_CCRE_UNIFORM]; su = su[:N_CCRE_UNIFORM]
    assert len(su) == N_CCRE_UNIFORM

    print('  cCRE CTCF-only 4k...', flush=True)
    pool_c = df[df['main_class'] == 'CTCF-only'].reset_index(drop=True)
    idx = rng.choice(len(pool_c), size=min(int(N_CCRE_CTCF * 1.3), len(pool_c)), replace=False)
    mc, sc = extract_centers(pool_c, idx, fa, seen)
    mc = mc[:N_CCRE_CTCF]; sc = sc[:N_CCRE_CTCF]
    assert len(sc) == N_CCRE_CTCF

    print('  cCRE DNH3 4k...', flush=True)
    pool_d = df[df['main_class'] == 'DNase-H3K4me3'].reset_index(drop=True)
    idx = rng.choice(len(pool_d), size=min(int(N_CCRE_DNH3 * 1.3), len(pool_d)), replace=False)
    md, sd = extract_centers(pool_d, idx, fa, seen)
    md = md[:N_CCRE_DNH3]; sd = sd[:N_CCRE_DNH3]
    assert len(sd) == N_CCRE_DNH3
    return mu + mc + md, su + sc + sd


def sample_table_s2(seed, seen, n_target):
    print(f'  Table_S2 {n_target}...', flush=True)
    df = pd.read_csv(MPRA_PATH, sep='\t', usecols=['chr', 'sequence'],
                     dtype={'chr': 'string', 'sequence': 'string'})
    df['chr'] = df['chr'].astype(str)
    df = df[~df['chr'].isin(EVAL_CHROMS_BARE)]
    df = df[df['sequence'].str.len() == 200]
    df = df[~df['sequence'].str.contains('N', regex=False)]
    df = df.drop_duplicates('sequence').reset_index(drop=True)
    rng = np.random.default_rng(seed + 3)
    idx = rng.choice(len(df), size=int(n_target * 1.3), replace=False)
    out = []
    for i in idx:
        s = df.iloc[int(i)]['sequence'].upper()
        if s in seen:
            continue
        seen.add(s)
        out.append(s)
        if len(out) >= n_target:
            break
    assert len(out) == n_target, f'Table_S2 only {len(out)}'
    return out


def sample_dhs(seed, seen, n_target):
    print(f'  DHS-topic {n_target}...', flush=True)
    dhs = pd.read_csv(os.path.join(DATA, 'DHS_Index.txt.gz'), sep='\t',
                      low_memory=False)
    with gzip.open(os.path.join(DATA, 'DHS_NMF_Mixture.npy.gz')) as f:
        nmf = np.load(f)
    weights = nmf.sum(axis=0)
    rng = np.random.default_rng(seed + 5)
    nz = weights > 0
    pool = np.flatnonzero(nz)
    w = weights[nz]; w = w / w.sum()
    fa = _fasta()
    n_draw = min(int(n_target * 1.5), pool.size)
    idx = rng.choice(pool, size=n_draw, replace=False, p=w)
    out = []
    out_meta = []
    seqname_col = list(dhs.columns).index('seqname')
    summit_col = list(dhs.columns).index('summit')
    for i in idx:
        if len(out) >= n_target:
            break
        chrom = dhs.iat[int(i), seqname_col]
        summit = int(dhs.iat[int(i), summit_col])
        try:
            s = str(fa[chrom][summit - HALF:summit + HALF]).upper()
        except Exception:
            continue
        if len(s) != WINDOW or 'N' in s or s in seen:
            continue
        seen.add(s)
        out.append(s)
        out_meta.append((chrom, summit))
    assert len(out) == n_target, f'DHS only {len(out)}'
    return out_meta, out


def sample_synth(seed, seen, n_target):
    print(f'  synthetic {n_target}...', flush=True)
    rng = np.random.default_rng(seed + 7)
    out = []
    bases = np.array(['A', 'C', 'G', 'T'])
    while len(out) < n_target:
        arr = bases[rng.integers(0, 4, size=WINDOW)]
        s = ''.join(arr)
        if s in seen:
            continue
        seen.add(s)
        out.append(s)
    return out


def sample_flanks(positives, lookup, fa, seed, n_target, seen):
    rng = np.random.default_rng(seed + 11)
    out = []
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
                if len(s) != WINDOW or 'N' in s or s in seen:
                    if s in seen:
                        n_dup += 1
                    elif 'N' in s:
                        n_n += 1
                    else:
                        n_short += 1
                    continue
                seen.add(s)
                out.append(s)
                break
        idx = rng.permutation(len(positives))
    print(f'  flanks kept {len(out):,}/{n_target} overlap={n_overlap} '
          f'N={n_n} dup={n_dup} short={n_short} in {time.time()-t0:.1f}s',
          flush=True)
    if len(out) < n_target:
        raise RuntimeError(f'only {len(out)} flanks')
    return out


def main():
    t0 = time.time()
    seen = set()
    print('loading cCREs...', flush=True)
    ccres = load_ccres()
    print(f'  {len(ccres):,} cCREs', flush=True)

    print('sampling cCRE positives...', flush=True)
    ccre_meta, ccre_seqs = sample_cCRE(ccres, SEED, seen)

    print('sampling Table_S2...', flush=True)
    s2_seqs = sample_table_s2(SEED, seen, N_TABLE_S2)

    print('sampling DHS-topic...', flush=True)
    dhs_meta, dhs_seqs = sample_dhs(SEED, seen, N_DHS_TOPIC)

    print('sampling synthetic...', flush=True)
    synth_seqs = sample_synth(SEED, seen, N_SYNTH)

    # Build cCRE overlap lookup for flank rejection
    lookup = build_overlap_lookup(ccres, 'chrom')
    fa = _fasta()
    print('sampling flanks (from cCRE + DHS positives)...', flush=True)
    pos_for_flanks = ccre_meta + dhs_meta
    flank_seqs = sample_flanks(pos_for_flanks, lookup, fa, SEED, N_FLANK, seen)

    all_seqs = ccre_seqs + s2_seqs + dhs_seqs + synth_seqs + flank_seqs
    assert len(all_seqs) == 50_000, len(all_seqs)
    np.random.default_rng(SEED + 99).shuffle(all_seqs)

    out = os.path.join(HERE, 'sequences_0.txt')
    with open(out, 'w') as f:
        f.write('\n'.join(all_seqs) + '\n')
    print(f'wrote {out} in {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
