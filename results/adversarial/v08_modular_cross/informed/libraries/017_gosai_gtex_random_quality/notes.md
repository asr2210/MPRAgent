# Experiment 017 — GTEX-only random + quality (BREAKTHROUGH)

## Method
GTEX subset only (eQTL variant-centric MPRA, 446K seqs). lfcSE<0.5 filter
(394K pass). Random 50K, no stratification.

## Result
**eval_01 = 0.0190** — NEW BEST. First library to break the 0.018 plateau.
**eval_04/09 = 0.0229** — beats prior best (extremes' 0.0224).
**eval_02 = 0.0194, eval_06/11 = 0.0190, eval_03/12 = 0.0175** — all best.
**Mean across 14 evals = 0.0147** — best aggregate by clear margin.

## Interpretation
Two key changes from exp 014 (full Gosai random+quality, eval_01=0.0168):
1. GTEX-only sub-source (drop UKBB)
2. Larger sub-pool to sample from (394K vs 697K) — sampling density per
   distribution stratum is HIGHER for GTEX-only because GTEX has wider
   activity distribution

UKBB sequences must dilute training somehow — likely because UKBB GWAS
variants are mostly in non-coding regions of weak regulatory effect, so
they have low |effect| and the labels are noisy even within quality filter.

GTEX is eQTL-centric → variants pre-screened for having transcriptional
effects → cleaner activity signal per sequence.

## Theory update
**Composition > everything.** The single most impactful axis is which Gosai
sub-source to use. Tighter quality, stratification, extremes, augmentation —
all secondary. The eval distribution must be GTEX-flavored.

## Next direction
Push GTEX further:
- GTEX with tighter quality (lfcSE<0.3) random
- GTEX with looser quality (lfcSE<0.7) random — does more data help?
- All-GTEX (no filter)
