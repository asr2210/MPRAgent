# 003_dhs_stratified

## Design
50k sequences from Meuleman SynthSeqs train set (160k pre-extracted
200bp DHS sequences, 10k per NMF component, chr3–chrY). 3,125 per
component, sampled without replacement, seed=0.

This is the published `dhs_stratified` strategy: 16 NMF cell-context
programs contribute equally regardless of their original DHS-pool size.

## Result
- eval_01 = **0.4378** (vs 002 genome_random: 0.4992; Δ -0.061)
- eval_07 = 0.4521 (vs 0.5985; Δ -0.146) — big drop on eval_07
- eval_13 = 0.4672 (vs 0.6020; Δ -0.135)
- eval_08 = 0.1016 (vs 0.0916; Δ +0.010, basically unchanged)
- mean over 14 ≈ 0.439 (vs 0.503; Δ -0.064)

## What this means — surprising contradiction with published table
The published baselines suggest DHS strategies should DOMINATE
synthetic (dhs_topic 0.7232 vs synth_oracle 0.6840). Translating to
my pipeline I expected DHS-stratified to be at or above my genome-
random baseline. Instead, DHS-stratified is uniformly worse.

Most likely explanation: in this pipeline, the eval distribution is
broader than the DHS distribution. Training only on DHS regions
narrows the model's exposure to the wider genomic prior (k-mer
distribution, repeats, gene bodies, etc.). The published-table
ordering presumably reflects a different model architecture or
evaluation distribution.

The pipeline gives me an answer that disagrees with the prior table:
**broad real-DNA dominates accessibility-restricted real-DNA** in
this evaluation. I'll keep trusting my own measurements over the
prior table.

## Confounds to rule out
- E3 uses chr3–chrY (SynthSeqs); E2 uses chr1/17/19/22. The
  comparison is partly cross-chromosome. Doesn't seem big enough to
  flip the sign of the effect, but worth a clean test.
- Single seed; differences <0.03 are noisy. Δ -0.061 on eval_01 is
  well outside noise. Δ -0.146 on eval_07 is huge.

## Generalization argument (now in doubt)
Original argument: stratifying across NMF components forces equal
exposure to all cell-lineage regulatory programs, which should
generalize across cell types. But this argument assumes the eval set
also lives in DHS regions. If the eval set is broader, the model
sees too narrow a distribution.

This reframes the problem: generalization to unseen cell types may
not require focusing on regulatory programs — it requires giving the
model the broadest possible exposure to the substrate (real DNA in
its natural distribution), so the model learns universal sequence
features rather than program-specific ones.
