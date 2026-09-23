# 02 — snRNA-seq quality control

Recorded 2026-09-23. Scripts: `src/02_qc_metrics_rna.py` (metrics and figure, no filtering),
`src/03_qc_filter_rna.py` (thresholds and doublet flags). Tables: `results/02_qc_rna_per_sample.tsv`,
`results/03_qc_rna_filter_summary.tsv`. Figure: `docs/figures/02_qc_rna_metrics.png`.

## What the metrics looked like before filtering

- Median UMI per nucleus ranged from 3,351 (Sample-52) to 8,344 (Sample-100); median genes detected from
  1,646 to 2,932. Depth differs by batch: the four Batch 1 samples are the deepest, Batch 3 controls the
  shallowest. This is the direct reason thresholds are set per sample.
- UMI and gene counts are bimodal in most samples, with a main mode near 6,000 UMI and a second mode near
  60,000 UMI. Of the 3,616 nuclei above log1p(UMI) = 10.3, 97.5 % carry the authors' EX or INH label
  (excitatory or inhibitory neuron); median log1p(UMI) is 10.13 for EX vs 8.67 for ODC. The upper mode is
  neurons, not doublets. **No upper bound is applied to counts or genes**; a symmetric MAD rule would remove
  a cell class.
- Mitochondrial fraction is near zero in nuclei as expected (medians 0.02 % to 0.34 %), with a thin tail
  reaching 20 % to 80 % in a few hundred nuclei. Five samples (27, 82, 19, 33, 90) sit a step higher than
  the rest (median 0.14 % to 0.34 %, 95th percentile 1.7 % to 2.6 %), consistent with slightly more
  cytoplasmic carry-over. A MAD rule is meaningless on a distribution centred this close to zero, so a fixed
  ceiling is used instead.

## Rules applied

| Rule | Threshold | Scope |
|---|---|---|
| Low UMI | log1p(total_counts) < median − 5 MAD | per sample |
| Low genes | log1p(n_genes) < median − 5 MAD | per sample |
| High mitochondrial | pct_counts_mt > 5 % | fixed, all samples |
| Doublet | scrublet predicted_doublet, run separately per sample, random_state 0 | per sample |

## Result

59,832 of 61,770 nuclei kept (96.9 %). Per-sample removal 1.0 % to 6.9 %. The mitochondrial rule removed
the most nuclei in exactly the five samples identified above (82, 55, 50, 49, 34 nuclei), so two independent
views agree on which samples carry more cytoplasmic contamination.

Agreement with the authors' own filtering: of the 298 nuclei the authors dropped, this pipeline drops 180
(60 %) and would keep 118; of the 61,472 the authors kept, this pipeline drops 1,758, almost all by the
mitochondrial and doublet rules, which the authors did not apply to this matrix.

## Known limitation: doublet calls are probably incomplete

Scrublet flagged 0.1 % to 2.7 % of nuclei per sample. For 10x loadings that recover 3,000 to 5,000 nuclei,
the expected doublet rate is roughly 3 % to 4 %, and the four Batch 1 samples with similar loading were
called at 3, 2, 55 and 38 doublets, which is not consistent. The likely cause is scrublet's automatic
threshold, which needs a bimodal score distribution and is known to be unstable on nuclei.

Decision: do not replace the automatic threshold with a hand-picked one, which would be equally arbitrary.
Keep `doublet_score` in `obs`, and after clustering check explicitly for clusters that co-express markers of
two lineages (for example a cluster positive for both MBP and RBFOX3); such clusters are removed at that
stage with the reason recorded.

## Output

`data/processed/rna_qc.h5ad`: 59,832 nuclei × 58,721 genes, raw counts in `X` and duplicated in
`layers["counts"]`, QC metrics, all four flags and `doublet_score` in `obs`.
