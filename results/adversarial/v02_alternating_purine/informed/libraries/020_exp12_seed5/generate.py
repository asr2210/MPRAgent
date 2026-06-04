"""012_dhs_celltype_matched — DHS enriched for the 3 model cell types.

DHS_Index has a 'component' column with chromatin program labels matching
the 3 labeled cell types:
- K562 (erythroid leukemia) -> Myeloid / erythroid (186k DHS)
- HepG2 (liver/hepatocyte)  -> Digestive (144k DHS)
- SK-N-SH (neuroblastoma)   -> Neural (461k DHS)

K562/HepG2 correlations have been near zero across all 11 experiments
(SKNSH dominates). Hypothesis: cell-type-matched DHS sequences provide
direct training signal for each cell-type head, lifting K562 and HepG2 r.

Composition (50k):
- 8k DHS Neural (SK-N-SH anchor)
- 6k DHS Myeloid/erythroid (K562 anchor)
- 6k DHS Digestive (HepG2 anchor)
- 5k cCRE uniform (general regulatory grammar)
- 5k Table_S2 random non-eval (variant MPRA coverage)
- 20k paired flanks (negatives, anchored to cCRE + DHS positives)
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
SEED = 5
FLANK_MIN = 1500
FLANK_MAX = 3000
EVAL_CHROMS_BARE = {'7', '13', '19', '21', 'X'}
MPRA_PATH = '/data/users/arao/.private/mpra_exp/data/Table_S2__MPRA_dataset.txt'

N_DHS_NEURAL = 8_000
N_DHS_MYELOID = 6_000
N_DHS_DIGESTIVE = 6_000
N_CCRE_UNI = 5_000
N_S2 = 5_000
N_FLANK = 20_000
assert (N_DHS_NEURAL + N_DHS_MYELOID + N_DHS_DIGESTIVE
        + N_CCRE_UNI + N_S2 + N_FLANK) == 50_000


def load_ccres():
    df = pd.read_csv(os.path.join(DATA, 'GRCh38-cCREs.bed'), sep='\t',
                     header=None,
                     names=['chrom', 'start', 'end', 'rDHS', 'accession', 'classes'])
    df['main_class'] = df['classes'].str.split(',').str[0]
    df['center'] = ((df['start'] + df['end']) // 2).astype(int)
    return df


def build_overlap_lookup(df):
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


def extract_centers(df, indices, fa, seen):
    seqs, meta = [], []
    for i in indices:
        row = df.iloc[int(i)]
        ws, we = row['center'] - HALF, row['center'] + HALF
        try:
            s = str(fa[row['chrom']][ws:we]).upper()
        except Exception:
            continue
        if len(s) != WINDOW or 'N' in s or s in seen:
            continue
        seen.add(s); seqs.append(s); meta.append((row['chrom'], int(row['center'])))
    return meta, seqs


def sample_class(df, klass, n, seed, seen):
    pool = df if klass is None else df[df['main_class'] == klass]
    pool = pool.reset_index(drop=True)
    rng = np.random.default_rng(seed)
    fa = _fasta()
    idx = rng.choice(len(pool), size=min(int(n * 1.4), len(pool)), replace=False)
    meta, seqs = extract_centers(pool, idx, fa, seen)
    return meta[:n], seqs[:n]


def _load_dhs():
    return pd.read_csv(os.path.join(DATA, 'DHS_Index.txt.gz'), sep='\t',
                       low_memory=False,
                       usecols=['seqname', 'summit', 'component', 'mean_signal'])


def sample_dhs_component(dhs, comp, n, seed, seen):
    """Sample n DHS sequences from given component, weighted by mean_signal."""
    sub = dhs[dhs['component'] == comp].reset_index(drop=True)
    chs = sub['seqname'].to_numpy()
    sus = sub['summit'].to_numpy()
    w = sub['mean_signal'].to_numpy()
    w = w / w.sum()
    rng = np.random.default_rng(seed)
    fa = _fasta()
    n_draw = min(int(n * 1.5), len(sub))
    idx = rng.choice(len(sub), size=n_draw, replace=False, p=w)
    seqs, meta = [], []
    for i in idx:
        if len(seqs) >= n:
            break
        ch = str(chs[int(i)]); su = int(sus[int(i)])
        try:
            s = str(fa[ch][su - HALF:su + HALF]).upper()
        except Exception:
            continue
        if len(s) != WINDOW or 'N' in s or s in seen:
            continue
        seen.add(s); seqs.append(s); meta.append((ch, su))
    assert len(seqs) == n, f'{comp} only {len(seqs)}/{n}'
    return meta, seqs


def sample_table_s2(seed, seen, n):
    df = pd.read_csv(MPRA_PATH, sep='\t', usecols=['chr', 'sequence'],
                     dtype={'chr': 'string', 'sequence': 'string'})
    df['chr'] = df['chr'].astype(str)
    df = df[~df['chr'].isin(EVAL_CHROMS_BARE)]
    df = df[df['sequence'].str.len() == 200]
    df = df[~df['sequence'].str.contains('N', regex=False)]
    df = df.drop_duplicates('sequence').reset_index(drop=True)
    rng = np.random.default_rng(seed + 3)
    idx = rng.choice(len(df), size=int(n * 1.3), replace=False)
    out = []
    for i in idx:
        s = df.iloc[int(i)]['sequence'].upper()
        if s in seen:
            continue
        seen.add(s); out.append(s)
        if len(out) >= n:
            break
    assert len(out) == n
    return out


def sample_flanks(positives, lookup, fa, seed, n, seen):
    rng = np.random.default_rng(seed + 11)
    out = []
    idx = rng.permutation(len(positives))
    pass_no = 0; t0 = time.time()
    while len(out) < n and pass_no < 8:
        pass_no += 1
        for i in idx:
            if len(out) >= n:
                break
            chrom, center = positives[int(i)]
            for _ in range(4):
                sign = 1 if rng.random() < 0.5 else -1
                off = sign * int(rng.integers(FLANK_MIN, FLANK_MAX + 1))
                c2 = center + off
                ws, we = c2 - HALF, c2 + HALF
                if ws < 0 or overlaps(chrom, ws, we, lookup):
                    continue
                try:
                    s = str(fa[chrom][ws:we]).upper()
                except Exception:
                    continue
                if len(s) != WINDOW or 'N' in s or s in seen:
                    continue
                seen.add(s); out.append(s); break
        idx = rng.permutation(len(positives))
    print(f'  flanks: {len(out):,}/{n} ({time.time()-t0:.1f}s)', flush=True)
    if len(out) < n:
        raise RuntimeError(f'only {len(out)}/{n} flanks')
    return out


def main():
    t0 = time.time()
    seen = set()
    print('loading cCREs...', flush=True)
    ccres = load_ccres()
    print('loading DHS_Index...', flush=True)
    dhs = _load_dhs()

    print(f'DHS Neural {N_DHS_NEURAL}...', flush=True)
    m_n, s_n = sample_dhs_component(dhs, 'Neural', N_DHS_NEURAL, SEED + 21, seen)
    print(f'DHS Myeloid/erythroid {N_DHS_MYELOID}...', flush=True)
    m_m, s_m = sample_dhs_component(dhs, 'Myeloid / erythroid', N_DHS_MYELOID, SEED + 22, seen)
    print(f'DHS Digestive {N_DHS_DIGESTIVE}...', flush=True)
    m_d, s_d = sample_dhs_component(dhs, 'Digestive', N_DHS_DIGESTIVE, SEED + 23, seen)
    print(f'cCRE uniform {N_CCRE_UNI}...', flush=True)
    m_u, s_u = sample_class(ccres, None, N_CCRE_UNI, SEED, seen)
    print(f'Table_S2 {N_S2}...', flush=True)
    s_s2 = sample_table_s2(SEED, seen, N_S2)

    lookup = build_overlap_lookup(ccres)
    fa = _fasta()
    pos_for_flanks = m_n + m_m + m_d + m_u
    print(f'flanks {N_FLANK} from {len(pos_for_flanks)} positives...', flush=True)
    s_f = sample_flanks(pos_for_flanks, lookup, fa, SEED, N_FLANK, seen)

    all_seqs = s_n + s_m + s_d + s_u + s_s2 + s_f
    assert len(all_seqs) == 50_000, len(all_seqs)
    np.random.default_rng(SEED + 99).shuffle(all_seqs)

    out = os.path.join(HERE, 'sequences_0.txt')
    with open(out, 'w') as f:
        f.write('\n'.join(all_seqs) + '\n')
    print(f'wrote {out} in {time.time()-t0:.1f}s', flush=True)


if __name__ == '__main__':
    main()
