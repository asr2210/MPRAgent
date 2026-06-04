# Skill 03 — Using the Malinois pretrained CNN as a sequence-activity oracle

## What it is
Malinois is the Gosai et al. 2024 BassetBranched CNN, pretrained on 776k MPRA
sequences from K562, HepG2, and SK-N-SH. Inference gives a `(B, 3)` tensor of
predicted log2FC for those three cell types.

It is a **prior model**, not the evaluator. Using it to pre-filter or
pre-score training candidates is legal and very fast (~15k seqs/s on 1 GPU).

## Setup
- Repo: `/data/users/arao/.private/mpra_exp/boda2`
- Venv: `/data/users/arao/.private/mpra_exp/venv/bin/python` (Torch 2.7.1+cu128)
- Pretrained checkpoint: `/data/users/arao/.private/mpra_exp/data/malinois_artifacts/artifacts/torch_checkpoint.pt`
- Constants: `boda.common.constants.MPRA_UPSTREAM` (300bp), `MPRA_DOWNSTREAM` (300bp)
- DON'T `import boda` (pulls dmslogo which is broken). Just import the
  submodules you need:
  ```python
  sys.path.insert(0, '/data/users/arao/.private/mpra_exp/boda2')
  from boda.model.basset import BassetBranched
  from boda.common.constants import MPRA_UPSTREAM, MPRA_DOWNSTREAM
  ```

## Loading the model
```python
import torch
ART = "/data/users/arao/.private/mpra_exp/data/malinois_artifacts/artifacts"
ckpt = torch.load(ART + "/torch_checkpoint.pt", map_location="cpu", weights_only=False)
hp = vars(ckpt["model_hparams"])
model = BassetBranched(**hp)
model.load_state_dict(ckpt["model_state_dict"])
model.eval().cuda()
```

## Input format
- (B, 4, 600) float one-hot encoding
- Encoding: A=0, C=1, G=2, T=3 → channel index
- Each 200bp core sequence must be padded with **200bp of MPRA flank on each
  side** for total length 600:
  - left: `MPRA_UPSTREAM[-200:]`
  - right: `MPRA_DOWNSTREAM[:200]`
- Output: `model(x)` returns (B, 3) tensor — column 0 = K562, 1 = HepG2,
  2 = SK-N-SH predicted log2FC.

## Fast batched scoring (~15k seqs/s on 1 H100/GB10)
Vectorize the one-hot:
```python
# core_idx: (B, 200) int8 in {0,1,2,3}
up_idx = np.array([m[c] for c in MPRA_UPSTREAM[-200:]], dtype=np.int8)
dn_idx = np.array([m[c] for c in MPRA_DOWNSTREAM[:200]], dtype=np.int8)
full = np.empty((B, 600), dtype=np.int64)
full[:, :200] = up_idx[None, :]
full[:, 200:400] = core_idx
full[:, 400:600] = dn_idx[None, :]
x = torch.zeros(B, 4, 600, device='cuda')
x.scatter_(1, torch.from_numpy(full).cuda().unsqueeze(1), 1.0)
```
- Use batch size 1024 (GPU memory friendly, throughput ≈ 15k/s)

## Output value scale (on random GC=50 sequences, n=500k)
- K562: range [-1.2, 7.9], mean 0.70
- HepG2: range [-1.2, 7.1], mean 0.65
- SK-N-SH: range [-1.5, 9.7], mean 0.67

Most sequences cluster near 0 (no signal); long right tail of high-activity
sequences. Quantile-based binning gives equal-count strata across the wide
dynamic range.

## Selection strategies
1. **Activity-span (3D)**: bin each axis into 5 quantiles → 125 cells →
   sample ~400 per cell. Forces broad coverage of joint predicted-activity
   space.
2. **Cell-type-discriminating**: select sequences with high variance across
   the 3 cell-type predictions (top-k by `var(preds, axis=1)`).
3. **High-activity-only**: top-k by mean across cell types (saturates the
   evaluator if activity range is the bottleneck).
4. **Adversarial**: combine an ensemble of priors and pick sequences where
   priors disagree → most informative for training.

## Caveats
- Malinois may have systematic blind spots → don't trust absolute log2FC,
  only relative ranks/quantiles.
- Predictions are conditioned on the MPRA_UPSTREAM/DOWNSTREAM flanks; the
  evaluator may use different (or no) flanks. Filtering on Malinois ranks
  rather than absolute predictions is more robust.
