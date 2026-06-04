# Experiment 005 — Table_S2 Malinois MPRA (50k random, no eval chroms)

## Design
- 50k random Table_S2 sequences (798k → 640k non-eval-chrom 200bp pool)
- Excludes chr7, 13, 19, 21, X (eval chromosomes per eval_mapping.json)
- No flanks; Table_S2 sequences are already MPRA-measured (mix active/inactive)

## Result — eval_01 = 0.1351
- K562_r = -0.026 / SKNSH = 0.457

## Major finding — eval-set-specific behavior
Per-eval comparison vs prior best:
- **eval_07 (sei_chr7_13): 0.193** — NEW HIGH for me (prior 0.16)
- **eval_13 (genomic_chr7_13): 0.166** — improved (prior 0.144)
- **eval_10 (dhs_chr7_13): 0.155** — slight improvement
- eval_06 (ukbb_gtex_both_oracle): 0.106 — DROPPED (prior 0.159)
- eval_11 (ukbb_gtex_both): 0.106 — DROPPED (prior 0.159)

Counterintuitive: Table_S2 is mostly UKBB/GTEx and EVAL_06/11 are UKBB/GTEx
"both alleles" evals. Despite source-relevance, Table_S2 alone HURT these.
Likely because the model needs paired ref/alt context or specific
distributional features the eval requires.

## Theory v4 — source diversity covers eval diversity
**No single source is uniformly best.** Each pool's distributional bias
excels on some evals and fails on others:
- DHS-topic → universal mediocre
- cCRE asymm → strong on natural eval_06/11
- Table_S2 → strong on eval_07/10/13 (chr-aware)
- The harness uses 14 different evaluation distributions; a library that
  matches NONE of them perfectly but covers the union will win.

## Plan
Exp 006 = direct replication of blind 013 (15k uniform cCRE + 5k CTCF +
5k DNH3 + 25k paired flanks). I need a pipeline anchor at the best-known
recipe to validate I can reach 0.17+ before designing combinations.
