# E002 — DHS uniform random (real labels)

## Design
50k sequences sampled uniformly at random from the Meuleman 2020 DHS index
(~3.59M sites, hg38). 200bp window centered on DHS summit. Seed=1.

## Result
eval_01 = **0.3179** (baseline `dhs_random` oracle-labeled: 0.7089 — 0.39 gap!)
eval_08 = 0.0797 (baseline oracle: 0.6673)
Per cell type eval_01: K562=0.143, HepG2=0.194, SKNSH=0.617
Total time: 76s (39s prepare eval + 37s extraction/training overhead)

## Key takeaways
1. **Real-vs-oracle gap is huge even for DHS**: 0.32 vs 0.71. Confirms baselines
   in instructions.md are predominantly under ORACLE LABELS, not real MPRA.
   My pipeline scores systematically lower across all strategies.

2. **SKNSH dominates across libraries**: random → 0.638, DHS → 0.617. K562/HepG2
   scores are weak (~0.14-0.19) under both. SKNSH might be easier to predict
   from sequence alone, or the eval cell types include several SKNSH-related ones.

3. **DHS only marginally beats random under real labels**: 0.32 vs 0.24. The
   strong DHS signal in oracle baselines may largely come from oracle's clean
   labeling, not from DHS sequences being inherently more learnable.

4. **mpra_real baseline at 0.6026** dramatically outperforms DHS-random with
   real labels (0.32). This suggests sequences from REAL MPRA datasets carry
   much more learnable signal than randomly sampled DHS sites — likely because
   MPRA datasets pre-select for sequences with measurable activity.

## Implications for next experiments
- The "informed" baseline table is largely IRRELEVANT to my real-label pipeline
- I need to find sequences with HIGH measurable activity (not just accessible)
- Should test: top-accessibility DHS, TF-bound DHS, published MPRA sequences,
  cCREs, ENCODE-validated enhancers

## Next experiment hypothesis
HIGH-SIGNAL DHS sites (top accessibility, multi-sample) should produce stronger
MPRA signal than uniform random DHS → real-label score should jump significantly.
If this works, it'll be a clean experiment showing "active sequences matter"
under real labels.
