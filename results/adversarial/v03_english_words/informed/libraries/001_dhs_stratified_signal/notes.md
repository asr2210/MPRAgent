# 001 — DHS NMF×signal-stratified

## Plan
Sample 50k from Meuleman 2020 synthseqs (160k pre-extracted 200bp DHS sequences).
Stratify equally across 16 NMF components AND 4 total_signal quartiles per
component (≈781 per stratum). Goal: span both regulatory program space AND
activity dynamic range, expecting to beat dhs_topic (0.7232).

## Result
eval_01 = 0.3883. **Catastrophic underperformance** vs dhs_topic baseline.
Mean across 14 evals: ~0.388. SK-N-SH correlation is ~0.06 across all evals;
K562/HepG2 ~0.54.

## Why it failed
The `train_all_classifier_light.csv.gz` from meuleman.org is **NOT a representative
DHS sample**. It is curated for a *classification* task — selecting the most
component-specific elements (mean proportion = 0.87, median 0.94). 50% of
sequences are detected in ≤6 of 733 biosamples (very rare/specific elements).

So my library is **dominated by cell-type-specific elements that are NOT active in
K562/HepG2/SK-N-SH** — the model can't learn from these because there's no signal.
SK-N-SH is hit hardest (0.06) likely because few of the 16 NMF components map
strongly to SK-N-SH-like contexts.

Equal-per-component stratification compounds this by forcing equal representation
across components that may have low coverage in our three measured cell types.

## What this teaches
1. **The synthseqs "classifier_light" subset is wrong for general MPRA library
   training.** It's a discriminator-training set, not a representative sampler.
2. **Activity-stratified sampling can be ACTIVELY HARMFUL if "low signal" means
   "signal in this dataset's cell types not the ones we measure".**
3. **A library only teaches the model about regions where the labeling cell types
   actually have signal.** Cell-type-specific regions from OTHER tissues are dead
   weight (low or zero signal in K562/HepG2/SK-N-SH).
4. **For generalization to unseen cell types**: the library still needs the
   labeling cell types to *produce informative measurements*. Diverse-but-silent
   sequences teach nothing.

## Next
Need a representative DHS pool that includes constitutively active elements
(active across many biosamples) AND cell-type-specific ones. Options:
- Download the full DHS Index (3.6M elements) + extract sequences from hg38
- Use sequences from the synthseqs subset but BIASED toward `numsamples` high
  (constitutive elements) which should be active in many cell types including
  ours
- Mix with ENCODE cCREs or other broader annotations

Cheapest experiment 2: re-use the same 160k pool but sample BIASED toward high
numsamples (more biosamples → more likely active in K562/HepG2/SK-N-SH).
