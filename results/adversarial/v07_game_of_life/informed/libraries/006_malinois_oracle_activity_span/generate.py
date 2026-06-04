"""Experiment 006 — Malinois-oracle activity-spanning library.

Generate ~500k GC=50 random ACGT candidates, score each with the pretrained
Malinois CNN (predicts K562, HepG2, SK-N-SH log2FC), then select 50k that
uniformly span the 3D predicted activity space.

Hypothesis (theory v4): if the 0.40 ceiling reflects training-set ACTIVITY
range (not sequence-diversity), selecting sequences with diverse predicted
activities should break the ceiling. If composition is the ONLY lever, this
should score ~0.39 like every other GC-50 library.

Malinois is the Gosai et al. 2024 BassetBranched CNN — pretrained on 776k
MPRA sequences from K562/HepG2/SK-N-SH. The predicted activities are imperfect
but rank-correlated with real activities; selecting on them should enrich
for real high-/low-activity training examples without needing labels.

Generalization justification: a library with broad predicted activity in the
3 labeled cell types contains sequences that activate/repress through cell-
type-shared mechanisms (TFs are shared across cell types; only their
expression levels differ). Selecting for predicted activity range in our
labeling cells is a proxy for activity range in any held-out cell type.
"""
from __future__ import annotations
import os
import sys
import time
import numpy as np
import torch

THIS = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.dirname(os.path.dirname(THIS))
sys.path.insert(0, REPO)
sys.path.insert(0, "/data/users/arao/.private/mpra_exp/boda2")

from utils.seqlib import write_sequences, SEQLEN
from boda.model.basset import BassetBranched
from boda.common.constants import MPRA_UPSTREAM, MPRA_DOWNSTREAM

OUT = os.path.join(THIS, "sequences_0.txt")
SCORES_OUT = os.path.join(THIS, "scores.npz")
SEED = 1
N_CAND = 500_000
N_TARGET = 50_000
N_BINS = 5      # per cell-type axis → 5^3 = 125 cells
BATCH = 1024
ART_DIR = "/data/users/arao/.private/mpra_exp/data/malinois_artifacts/artifacts"


def load_malinois(device):
    ckpt = torch.load(os.path.join(ART_DIR, "torch_checkpoint.pt"),
                       map_location="cpu", weights_only=False)
    hp = vars(ckpt["model_hparams"])
    model = BassetBranched(**hp)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    return model.to(device)


def gen_gc50_indices(rng, n):
    """Generate n random 200bp sequences as int8 arrays (A=0,C=1,G=2,T=3). Shape (n, 200)."""
    return rng.integers(0, 4, size=(n, SEQLEN), dtype=np.int8)


def flank_indices():
    """Return (200,) int arrays for the left and right MPRA flanks."""
    m = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 0}  # N (rare) → A; flanks have no N anyway
    up = MPRA_UPSTREAM[-200:]
    dn = MPRA_DOWNSTREAM[:200]
    return np.array([m[c] for c in up], dtype=np.int8), np.array([m[c] for c in dn], dtype=np.int8)


_UP_IDX, _DN_IDX = None, None
def get_flanks():
    global _UP_IDX, _DN_IDX
    if _UP_IDX is None:
        _UP_IDX, _DN_IDX = flank_indices()
    return _UP_IDX, _DN_IDX


def onehot_batch_from_idx(core_idx, device):
    """One-hot encode core int indices padded with flanks. (B, 4, 600)."""
    up, dn = get_flanks()
    B = core_idx.shape[0]
    full = np.empty((B, 600), dtype=np.int64)
    full[:, :200] = up[None, :]
    full[:, 200:400] = core_idx
    full[:, 400:600] = dn[None, :]
    x = torch.zeros(B, 4, 600, device=device)
    full_t = torch.from_numpy(full).to(device)
    x.scatter_(1, full_t.unsqueeze(1), 1.0)
    return x


def score_all(model, all_idx, device):
    N = all_idx.shape[0]
    preds = np.zeros((N, 3), dtype=np.float32)
    t0 = time.time()
    for i in range(0, N, BATCH):
        batch_idx = all_idx[i:i + BATCH]
        x = onehot_batch_from_idx(batch_idx, device)
        with torch.no_grad():
            y = model(x)
        preds[i:i + BATCH] = y.cpu().numpy()
        if (i // BATCH) % 50 == 0:
            dt = time.time() - t0
            done = i + batch_idx.shape[0]
            rate = done / max(dt, 1e-6)
            eta = (N - done) / max(rate, 1e-6)
            print(f"  scored {done:,}/{N:,} | {rate:.0f}/s | ETA {eta:.0f}s", flush=True)
    return preds


def idx_to_seqs(idx_arr):
    """Convert (n, 200) int8 array to list of strings."""
    lookup = np.array(list("ACGT"))
    chars = lookup[idx_arr]  # (n, 200) str
    return ["".join(row) for row in chars]


def stratified_select(preds, rng, n_target, n_bins):
    """Bin sequences into n_bins**3 cells by per-axis quantile, then sample.

    For each axis, compute quantile boundaries from the full distribution so
    bins have equal counts. Then pick ~ n_target/n_bins^3 from each cell. If a
    cell is short, take all; if long, sample uniformly.
    """
    n = preds.shape[0]
    edges = np.quantile(preds, np.linspace(0, 1, n_bins + 1)[1:-1], axis=0)
    # edges shape: (n_bins-1, 3)
    bins = np.zeros((n, 3), dtype=np.int8)
    for k in range(3):
        bins[:, k] = np.digitize(preds[:, k], edges[:, k])  # 0..n_bins-1
    keys = bins[:, 0] * (n_bins * n_bins) + bins[:, 1] * n_bins + bins[:, 2]
    total_cells = n_bins ** 3
    target_per = n_target // total_cells
    leftover = n_target - target_per * total_cells

    selected = []
    cell_counts = np.bincount(keys, minlength=total_cells)
    print(f"  Cell count stats: min={cell_counts.min()}, median={int(np.median(cell_counts))}, max={cell_counts.max()}, n_empty={int((cell_counts==0).sum())}")
    cells_order = rng.permutation(total_cells)
    extras_per_cell = np.zeros(total_cells, dtype=int)
    for c in cells_order[:leftover]:
        extras_per_cell[c] = 1

    for c in range(total_cells):
        idxs = np.where(keys == c)[0]
        want = target_per + extras_per_cell[c]
        if len(idxs) == 0:
            continue
        if len(idxs) <= want:
            selected.extend(idxs.tolist())
        else:
            pick = rng.choice(idxs, size=want, replace=False)
            selected.extend(pick.tolist())

    selected = np.array(selected)
    # Top off if undersized due to empty cells: random pick from remaining
    remaining = np.setdiff1d(np.arange(n), selected, assume_unique=False)
    if len(selected) < n_target:
        deficit = n_target - len(selected)
        extra = rng.choice(remaining, size=deficit, replace=False)
        selected = np.concatenate([selected, extra])
    elif len(selected) > n_target:
        keep = rng.choice(len(selected), size=n_target, replace=False)
        selected = selected[keep]
    rng.shuffle(selected)
    return selected


def main():
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    print(f"Device: {device}")

    print("Loading Malinois...")
    model = load_malinois(device)

    # Generate candidates as int8 indices (compact, vectorizable)
    print(f"Generating {N_CAND:,} GC-50 random candidates...", flush=True)
    all_idx = gen_gc50_indices(rng, N_CAND)
    print(f"  Generated shape {all_idx.shape} in {time.time()-t0:.1f}s", flush=True)

    print(f"Scoring with Malinois (batch={BATCH})...", flush=True)
    preds = score_all(model, all_idx, device)
    print(f"Pred ranges: K562 [{preds[:,0].min():.2f}, {preds[:,0].max():.2f}],"
          f" HepG2 [{preds[:,1].min():.2f}, {preds[:,1].max():.2f}],"
          f" SKNSH [{preds[:,2].min():.2f}, {preds[:,2].max():.2f}]")
    print(f"  means: K562 {preds[:,0].mean():.3f}, HepG2 {preds[:,1].mean():.3f}, SKNSH {preds[:,2].mean():.3f}")

    print(f"Stratified selection into {N_BINS}^3={N_BINS**3} cells...")
    sel = stratified_select(preds, rng, N_TARGET, N_BINS)
    print(f"  Selected {len(sel):,} sequences")

    selected_seqs = idx_to_seqs(all_idx[sel])
    selected_preds = preds[sel]
    print(f"Selected pred ranges: K562 [{selected_preds[:,0].min():.2f}, {selected_preds[:,0].max():.2f}],"
          f" HepG2 [{selected_preds[:,1].min():.2f}, {selected_preds[:,1].max():.2f}],"
          f" SKNSH [{selected_preds[:,2].min():.2f}, {selected_preds[:,2].max():.2f}]")
    print(f"Selected stds: K562 {selected_preds[:,0].std():.3f}, HepG2 {selected_preds[:,1].std():.3f}, SKNSH {selected_preds[:,2].std():.3f}")

    # GC stats
    gcs = np.array([sum(1 for c in s if c in "GC") / SEQLEN for s in selected_seqs])
    print(f"Selected GC: mean={gcs.mean():.3f}, std={gcs.std():.3f}")

    np.savez(SCORES_OUT,
              preds_all=preds.astype(np.float16),
              selected_idx=sel.astype(np.int32),
              selected_preds=selected_preds.astype(np.float16))
    write_sequences(selected_seqs, OUT)
    print(f"Done in {time.time()-t0:.1f}s")


if __name__ == "__main__":
    main()
