# Skill: cell-type-matched DHS enrichment

## What
The DHS_Index (`data/DHS_Index.txt.gz`) has a `component` column labeling each DHS by its chromatin program. For a model trained on K562/HepG2/SK-N-SH, sampling DHS sequences from matching components directly lifts the per-cell-type correlation r.

## Mapping
- K562 (erythroid leukemia) ↔ `Myeloid / erythroid` (186k DHS)
- HepG2 (liver) ↔ `Digestive` (144k DHS)
- SK-N-SH (neuroblastoma) ↔ `Neural` (461k DHS)

## How
```python
dhs = pd.read_csv('data/DHS_Index.txt.gz', sep='\t',
                  usecols=['seqname','summit','component','mean_signal'])
sub = dhs[dhs['component'] == 'Myeloid / erythroid']
# Weight by mean_signal (DHS activity strength). Sample, extract ±100bp around summit.
```

GOTCHA: `pd.read_csv(usecols=...)` does NOT preserve passed-in column order; it follows file order. Don't index by position — use column names or `to_numpy()`.

## Effect
+0.03 to +0.05 K562 r on UKBB/GTEx variant MPRA evals (eval_06/11, eval_03/12).
+0.011 mean_r(14), +0.015 eval_01.

## When to use
Whenever the target eval set evaluates cell-type-specific activity correlation. Cell-type-matched DHS gives the model direct per-head training signal that single-class cCREs and untargeted DHS-topic do not.

## Verified in
- libraries/012_dhs_celltype_matched (mean_r=0.162, eval_01=0.174)
