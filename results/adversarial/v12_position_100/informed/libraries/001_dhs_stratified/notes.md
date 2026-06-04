# 001 — DHS NMF-stratified

## What I built
50,000 DHS sequences from Meuleman et al. 2020 curated training subset (160k sequences, 10k per of 16 NMF components, chr3-chrY hg38). I sampled exactly 3125 from each NMF component → perfectly stratified across regulatory programs.

## What I predicted
This should be very close to the `dhs_stratified` strategy in instructions.md (eval_01 = 0.7055 there). At minimum, should clearly beat random_uniform (0.0711 in strategies.md).

## What happened
- eval_01 = 0.0739 — barely above random_uniform (0.0711), BELOW dirichlet_composition (0.0768).
- mean across evals ≈ 0.091. Best evals: eval_07 (0.1454), eval_13 (0.1408), eval_10 (0.1286). Worst: eval_08 (0.0627).
- 37s runtime.

## Interpretation
The huge gap between instructions.md (dhs_topic = 0.72 at 50k) and what I got here (0.07) means the two harnesses are very different. Instructions.md baselines almost certainly use oracle/teacher-model labels (which give signal regardless of library); THIS harness uses real MPRA measurements where signal-to-noise at 50k is brutal.

In the strategies.md regime, biological grounding gives at best a ~0.01-0.02 absolute boost over random — much smaller than the proportional gain in the oracle regime. The key bottleneck appears to be the noise floor of real measurements, not the choice of sequences per se.

That said: eval_07, eval_10, eval_13 systematically score higher for biological sequences (0.13-0.15) than for random (0.13ish in strategies.md too). These three evals appear more biology-sensitive. Worth tracking.

## What this updates in my theory
- H1 (motifs are core): not falsified, but the absolute gain from biological sequences over random is small. Need to figure out WHY.
- New hypothesis (H5): At 50k with real labels, the bottleneck is library-level diversity in ACTIVITY (dynamic range), not just sequence diversity. Random sequences with extreme compositional spread (dirichlet) might give the model more contrast to fit. Biological DHS sequences are functionally narrow (all open chromatin → all moderately active).
- This predicts: combining biologically-grounded sequences with extreme/contrast sequences should help more than either alone.

## What to try next
Several candidate experiments:
1. Half DHS + half random: does adding random contrast help?
2. SEI chromatin states: regions with very different expected activity (active enhancer vs heterochromatin) → would test the dynamic range hypothesis directly.
3. Strong promoters (high signal_total subset of DHS, top decile): bias toward HIGH-activity, large effect sizes.
4. Bulk download full DHS index, sample weighted by total_signal — get the full diversity.

Most informative seems to be (2): explicit dynamic range via chromatin state diversity.
