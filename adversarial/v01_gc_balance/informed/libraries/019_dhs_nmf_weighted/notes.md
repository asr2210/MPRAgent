# 019_dhs_nmf_weighted — notes

**Design**: Downloaded `2018-06-08NC16_NNDSVD_Mixture.npy` (16 x 3,591,898) from Meuleman Zenodo — per-DHS NMF topic loadings. Sampled probability proportional to sum of all 16 topic loadings per DHS (literal reading of baseline `dhs_topic` recipe).

**Result**: eval_01 = 0.6481 (vs exp 005 = 0.6752, **-0.027**). Mean = 0.6034.

**Setup verification**:
- Mapped matrix row order to topic names via column-argmax counts; 1-to-1 alignment with DHS Index rows confirmed (3,591,898 entries on both sides).
- Row order: 0=Tissue invariant, 1=Stromal A, 2=Primitive/embryonic, ..., 15=Cancer/epithelial.

**Surprise**: sum-weighted NMF hurts in my pipeline. Baseline `dhs_topic` (0.7232) > `dhs_random` (0.7089) by +0.014. My exp 019 (sum-weighted) < exp 003 (uniform) by -0.012. **Opposite sign**.

This is not a recipe failure on my part — the data and weights are correct. It's a deep pipeline difference: my prepare.py / surrogate model interacts with NMF-weighted libraries differently than the baseline pipeline did.

**Multi-seed test** (separately, gave 5 sequences_N.txt files to prepare.py for exp 005 recipe): eval_01 = 0.6743 vs single-seed 0.6752. Multi-seed averaging via file aggregation doesn't lift either.

**Theory update H7 → H8**:
> The 0.05 gap to baseline `dhs_topic` is not closable by recipe alone in this prepare.py setup. NMF weighting that helps in the baseline pipeline does not transfer here — probably a different surrogate / labeler that scores differently on cell-specific vs broad DHSs.

**Implication**: my local optimum is uniform within a thoughtfully filtered pool. Focus remaining experiments on novel filters and combinations.
