---
slug: aucell-ucell
method: "AUCell -> UCell (gene signature scoring)"
omicverse_binding: "omicverse.single.aucell (vendored AUCell/ctxcore recovery curves); UCell itself is not exposed"
original:
  paper: "Aibar et al., Nature Methods 2017, doi:10.1038/nmeth.4463 (AUCell introduced in SCENIC)"
  repo: "https://github.com/aertslab/AUCell"
  commit: "91753b327a39dc7a4bbed46408ec2271c485f2f0"
optimized:
  paper: "Andreatta & Carmona, Comput Struct Biotechnol J 2021, doi:10.1016/j.csbj.2021.06.043"
  repo: "https://github.com/carmonalab/UCell"
  commit: "f778e091f19e28d3d05bde134663fd6b2f938267"
optimization_type: algorithmic
mechanism: "Replace the integrated area under each cell's rank-recovery curve with the algebraically identical normalized Mann-Whitney U statistic, computed as one clipped rank sum per gene set for all cells in a chunk."
reported_speedup: "~3x faster than AUCell and up to 10x lower peak memory; 100,000 cells need 5.5 GB peak in UCell versus >64 GB in AUCell (UCell paper Fig. 1E-F, Hao et al. T-cell dataset; hardware not reported)"
equivalence: approximate
optimization_families: [O1, O6, O7, O10]
sources:
  - "https://doi.org/10.1016/j.csbj.2021.06.043"
  - "https://pmc.ncbi.nlm.nih.gov/articles/PMC8271111/"
  - "https://github.com/carmonalab/UCell/blob/f778e091f19e28d3d05bde134663fd6b2f938267/R/HelperFunctions.R"
  - "https://github.com/carmonalab/UCell/blob/f778e091f19e28d3d05bde134663fd6b2f938267/R/StoreRankings_UCell.R"
  - "https://github.com/aertslab/AUCell/blob/91753b327a39dc7a4bbed46408ec2271c485f2f0/R/02_calcAUC.R"
  - "https://github.com/aertslab/AUCell/blob/91753b327a39dc7a4bbed46408ec2271c485f2f0/R/01_buildRankings.R"
verified: true
---

## Baseline bottleneck

AUCell ranks every gene in every cell (`DelayedMatrixStats::colRanks(..., ties.method="random", decreasing=TRUE)`, `R/01_buildRankings.R:191-195`) and then, for each gene set and each cell, integrates the rank-recovery curve up to `aucMaxRank` (default 5% of genes). The inner loop is a per-cell R closure over a column of the gene-set rank submatrix:

```r
# AUCell: R/02_calcAUC.R:275-282
.auc <- function(oneRanking, aucThreshold, maxAUC)
{
  x <- unlist(oneRanking)
  x <- sort(x[x<aucThreshold])     # per-cell sort
  y <- seq_along(x)
  sum(diff(c(x, aucThreshold)) * y)/maxAUC
}
```

Costs: (i) an `apply`-level R call plus a sort and two vector allocations per (cell, gene set); (ii) a per-cell ranking matrix retained for the whole dataset, whose footprint dominates peak memory (the paper measures >64 GB for 100k cells); (iii) rankings are recomputed on every call unless the user explicitly stores them.

## Optimization trajectory

1. **Drop the curve, keep a sufficient statistic** — the area under a monotone recovery curve up to `T` equals `k*T - sum(ranks)`, so the whole statistic is a single sum of the gene set's ranks (analytic reformulation). No sort, no per-threshold accumulation.
2. **Clip and sum once per gene set, vectorized over cells** — the rank submatrix rows for the signature are clipped at `maxRank`, then one `apply(..., 2, sum)` produces the score for every cell in the chunk:

```r
# UCell: R/HelperFunctions.R:32-45
if (length(present_idx) > 0) {
  rank_sub <- ranks_matrix[present_idx, -1, drop = FALSE]
  if (sparse==TRUE){ rank_sub[rank_sub==0] <- maxRank }
  insig <- rank_sub >= maxRank
  rank_sub[insig] <- maxRank
  rank_sum = rank_sum + apply(rank_sub, 2, sum)
}
rank_sum_min <- len_sig*(len_sig + 1)/2
ucell_score <- 1 - (rank_sum-rank_sum_min)/(len_sig*maxRank - rank_sum_min)
```

3. **Bound memory by chunking cells** — `chunk.size=100` cells are densified at a time and scored, so peak memory is `100 x n_genes` instead of `n_cells x n_genes`; chunks are distributed with `BiocParallel::bplapply` (stream/chunk,  parallelism).
4. **Make ranking reusable** — `StoreRankings_UCell` stores one sparse rank matrix (ranks `>= maxRank` collapsed to 0) that `ScoreSignatures_UCell(precalc.ranks=...)` reuses for any number of signatures (reuse target-invariant state). This is the step that removes repeated ranking when scanning many signatures.
5. **Ranks computed by `data.table::frankv`** on a column-wise `as.data.table` copy, instead of per-column R sorting.

Relevant locations: `R/HelperFunctions.R:18-48` (`u_stat`), `R/HelperFunctions.R:124-193` (`calculate_Uscore`, chunking + `bplapply`), `R/HelperFunctions.R:325-333` (`data_to_ranks_data_table`), `R/StoreRankings_UCell.R:1-60`.

## Evidence

- "Compared to AUCell, UCell is about three times faster (Fig. 1E) and uses significantly less memory (Fig. 1F). For example, AUCell requires over 64 GB of RAM to process 100,000 single-cells, while UCell uses only 5.5 GB of peak memory (Fig. 1F)" — https://pmc.ncbi.nlm.nih.gov/articles/PMC8271111/
- "UCell is about three times faster and requires up to ten times less memory on large (>10^4) single-cell datasets" (Fig. 1 caption, same URL).
- Abstract: scores "demand less computing time and memory than other available methods, enabling the processing of large datasets in a few minutes even on machines with limited computing power" — https://doi.org/10.1016/j.csbj.2021.06.043
- Hardware, cell/gene dimensions beyond "datasets of different sizes" and the specific signature count are **not reported** in the paper text; timings are in Fig. 1E without a machine specification.

## Equivalence and caveats

The two scores are not merely correlated: for the same ranking, the same gene set and `aucMaxRank == maxRank`, AUCell's `normAUC=TRUE` value equals UCell's score identically, because

`AUC = (k*T - sum(rank_i)) / (k*T - k(k+1)/2)` and `UCell = 1 - (sum(rank_i) - k(k+1)/2) / (k*T - k(k+1)/2)`.

So the rewrite is exact as a statistic; the entry is marked `approximate` because the shipped defaults differ on three points that change numbers and even cell ordering:

- **Threshold**: AUCell `aucMaxRank = ceiling(0.05 * n_genes)`; UCell `maxRank = 1500`. These are different truncations of the same curve.
- **Ties**: AUCell defaults to `ties.method="random"` while building rankings; UCell defaults to `"average"`. On sparse count matrices, zero ties dominate, so rankings — and scores — differ beyond Monte Carlo noise.
- **Missing genes**: UCell imputes absent signature genes at `maxRank` (`missing_genes="impute"`, default) or drops them (`"skip"`); AUCell drops them and normalizes by the number of *present* genes.
- UCell adds `w_neg` (positive/negative gene subsets) with no AUCell equivalent.

Downstream cell clustering/annotation must be re-validated whenever defaults are not explicitly matched; matching `maxRank`/`aucMaxRank`, tie method and missing-gene policy makes the mapping exact.

## Reuse plan for Bio-Opt

- **Test first **: check whether the scoring statistic is a sufficient statistic of the ranked vector (a sum, count or quantile) rather than an integral over a curve; if yes, replace the integration with the closed form and verify identity on a small matrix before timing.
- **Then ()**: bound memory by chunking cells and cache the one expensive shared object (the rank matrix) keyed on its full dependency set (matrix hash, tie method, `maxRank`), reusing it across all signatures.
- Acceptance criteria: cell-level score equality within `1e-12` against a matched-configuration AUCell run on the same ranking; downstream Spearman/Kendall on cell ordering unchanged; peak RSS measured as a separate metric.
