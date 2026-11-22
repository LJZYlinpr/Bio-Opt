---
slug: spark-sparkx
method: "SPARK -> SPARK-X (spatially expressed gene detection)"
omicverse_binding: "not bound in the inspected omicverse 2.4.3 tree: omicverse.space.svg exposes mode={'prost','moran','spateo','somde','spatialde','pearson_residuals'} and no SPARK/SPARK-X reference exists in omicverse/space, omicverse/utils or the bundled guides; verified by repository-wide grep (binding status: unverified for older releases)"
original:
  paper: "Sun S, Zhu J, Zhou X. Statistical analysis of spatial expression patterns for spatially resolved transcriptomic studies. Nat Methods 18:193-200 (2021), doi:10.1038/s41592-019-0701-7"
  repo: "https://github.com/xzhoulab/SPARK"
  commit: "a8b4bf27b804604dfda53da42992f100b8e4e727"
optimized:
  paper: "Zhu J, Sun S, Zhou X. SPARK-X: non-parametric modeling enables scalable and robust detection of spatial expression patterns for large spatial transcriptomic studies. Genome Biol 22:184 (2021), doi:10.1186/s13059-021-02404-0"
  repo: "https://github.com/xzhoulab/SPARK"
  commit: "a8b4bf27b804604dfda53da42992f100b8e4e727"
optimization_type: algorithmic
mechanism: "Replace per-gene penalized quasi-likelihood fitting of a Poisson mixed model (iterative AI-REML over an n-by-n covariance) with a closed-form distance-covariance test whose d-by-d eigenproblem and kernel spectrum are computed once and shared by all genes."
reported_speedup: "orders of magnitude: 1 min versus 770 min (SpatialDE) and 1125 min (SPARK-G) for 10,000 genes on 20,000 locations, single thread of an Intel Xeon E5-2683 2.00 GHz (Fig. 1B); 0.32 GB versus 21.2/55.3 GB RAM (Fig. 1C); 3 min versus 56 h (SPARK-G) and 47 h (SpatialDE) on mouse cerebellum Slide-seq whole data; 2 min versus 13 h/8 h on Slide-seqV2; 3 min and 0.42 GB on 177,455 HDST spots, where it was the only applicable method"
equivalence: changed-output
optimization_families: [O1, O6, O10, O11]
sources:
  - "https://doi.org/10.1038/s41592-019-0701-7"
  - "https://doi.org/10.1186/s13059-021-02404-0"
  - "https://pmc.ncbi.nlm.nih.gov/articles/PMC8218388/"
  - "https://github.com/xzhoulab/SPARK/blob/a8b4bf27b804604dfda53da42992f100b8e4e727/R/sparkx.R"
  - "https://github.com/xzhoulab/SPARK/blob/a8b4bf27b804604dfda53da42992f100b8e4e727/R/vcestimating.R"
  - "https://github.com/xzhoulab/SPARK/blob/a8b4bf27b804604dfda53da42992f100b8e4e727/src/spark_float.cpp"
  - "https://zenodo.org/record/4903349"
verified: true
---

## Baseline bottleneck

SPARK treats each gene as a Poisson (or Gaussian) mixed model with a spatially correlated random effect: `log E[y] = offset + X alpha + b`, `b ~ N(0, tau * K)`, where `K` is one of ten spatial kernels: five Gaussian and five periodic/cosine kernels at different length scales (`R/vctesting.R:30-52`). For every gene and every kernel it estimates the variance component `tau` by penalized quasi-likelihood: a Fisher-scoring/AI-REML loop (`spark.ai`) nested inside an outer loop that repeatedly re-fixes variance components (`while(any(fixtau.new != fixtau.old & count_num<maxiter))`, `maxiter=500`), with the heavy linear algebra in C++ (`src/spark_float.cpp`).

Work per gene therefore scales with (i) the number of iterations to convergence, (ii) an `n x n` dense distance/kernel matrix (`ED <- as.matrix(dist(object@location[,1:2]))`, `R/vctesting.R:32`) that must be factorised or multiplied, and (iii) the ten kernels at fixed length-scale quantiles from `ComputeGaussianPL`. This is why the SPARK-X paper can state that SPARK was not applied above 3,000 locations in its benchmark, and that SPARK "fails to converge for all genes" on sparse data.

## Optimization trajectory

1. **Replace the model by a covariance test .** SPARK-X no longer estimates a variance component. It tests whether expression `y` covaries with the spatial kernel: the statistic is a normalized `y^T H y` (a distance-covariance / KRV-style quantity) with `H` built from the centered, covariate-adjusted coordinates. Under the null the statistic is a weighted sum of `chi^2_1` variables, so no likelihood, no link function and no iteration are needed (paper, Methods; `R/sparkx.R`).
2. **Reduce the statistic to a d-by-d eigenproblem ().** Because the kernel is a linear/projection kernel in the coordinate space, `H` is never formed. The code centers the coordinates (`Xinfomat`), computes `loc_inv = solve(crossprod(Xinfomat, Xinfomat))` once, gets the kernel spectrum `Klam` from `eigen(crossprod(Xinfomat, Xinfomat %*% loc_inv), only.values=TRUE)`, and then evaluates all genes with one matrix product `EHL = counts %*% Xinfomat` (`R/sparkx.R:127-141`). The `d x d` eigen-decomposition is shared by every gene and every kernel; memory and time are linear in the number of locations and genes, not quadratic.
3. **Closed-form p-values.** For each gene the null distribution is the weighted sum `lambda_G * lambda_K`, and the p-value comes from Davies' exact method, with a Liu moment-matching fallback when Davies returns a non-positive value (`sparkx_pval`, `R/sparkx.R:235-249`). This replaces per-gene iterative fitting with one small numerical integration.
4. **Vectorise and thread over genes .** The statistic is a single `apply` over `EHL`, and the per-gene p-value loop runs under `mclapply` (`R/sparkx.R:138-142`); sparse count matrices are accepted directly and zero cells/genes are dropped once up front.
5. **Keep the multi-kernel combination, enlarge the kernel set.** Eleven kernels (projection, 5 Gaussian, 5 cosine) are tested and their p-values combined by the Cauchy rule (ACAT) with BY adjustment (`R/sparkx.R:86-93`); the original also combined ten kernels with ACAT, so this step adds pattern coverage without adding an extra fitting loop.

```r
# R/sparkx.R:127-141 (sparkx.sk, no-covariate branch)

Xinfomat <- apply(infomat, 2, scale, scale = FALSE)
loc_inv  <- solve(crossprod(Xinfomat, Xinfomat))
kmat_first <- Xinfomat %*% loc_inv
LocDim   <- ncol(infomat)
Klam     <- eigen(crossprod(Xinfomat, kmat_first), only.values = TRUE)$values
EHL      <- counts %*% Xinfomat
numCell  <- nrow(Xinfomat)
adjust_nominator <- as.vector(sp_sums_Rcpp(counts^2, TRUE))
vec_stat <- apply(EHL, 1, function(x) x %*% loc_inv %*% as.matrix(x)) * numCell / adjust_nominator
vec_ybar <- as.vector(sp_means_Rcpp(counts, TRUE))
vec_ylam <- unlist(mclapply(1:nrow(counts), function(x)
                1 - numCell * vec_ybar[x]^2 / adjust_nominator[x], mc.cores = mc_cores))
vec_daviesp <- unlist(mclapply(1:nrow(counts), function(x)
                sparkx_pval(x, vec_ylam, Klam, vec_stat), mc.cores = mc_cores))
```

The loop it replaces is the PQL recursion in `R/vcestimating.R:204-216`:

```r
fixtau.old <- rep(0, num_vc + 1)
model1 <- spark.ai(model0, num_vc, maxiter = maxiter, tol = tol, verbose = verbose)
fixtau.new <- 1 * (model1$theta < 1.01 * tol)
count_num <- 1
while (any(fixtau.new != fixtau.old & count_num < maxiter)) {
    count_num <- count_num + 1
    fixtau.old <- fixtau.new
    model1 <- spark.ai(model0, num_vc, fixtau = fixtau.old, maxiter = maxiter, tol = tol, verbose = verbose)
    fixtau.new <- 1 * (model1$theta < 1.01 * tol)
}
```

## Evidence

All timings from the SPARK-X paper, single thread of an Intel Xeon E5-2683 2.00 GHz, 10,000 genes:

- Simulation, 20,000 locations: "it takes SPARK-G and SpatialDE 1125 and 770 min ... In contrast, it only takes SPARK-X 1 min to analyze the same data" and "SPARK-G and SpatialDE require 21.2 and 55.3 gigabytes (GB) of physical RAM memory, respectively, SPARK-X only requires 0.32 GB" (Results, Fig. 1B-C).
- Mouse cerebellum Slide-seq: "it took SPARK-X 3 min to analyze the whole data, while it took 56 h for SPARK-G and 47 h for SpatialDE" (Fig. 2C, Table S1).
- Slide-seqV2: "it took SPARK-X 2 min ... while it took 13 h for SPARK-G and 8 h for SpatialDE".
- HDST olfactory bulb, 19,913 genes x 177,455 spots: "SPARK-X requires 0.42 GB memory and 3 min of computing time and is the only SE method applicable to the data."
- Scaling statement: "with moderate computation resource, SPARK-X is the only method applicable to data with sample size exceeding around 30,000."

## Equivalence and caveats

`changed-output`: SPARK-X does not fit the same generative model, so it does not reproduce SPARK's likelihood-ratio statistics, variance components, or gene lists. On human heart Visium the paper reports 1,536 SE genes for SPARK-X versus 644 for SPARK (533 overlapping); on ovarian cancer, 651 versus 579 (474 overlapping). The tests are calibrated and at least as powerful in the paper's benchmarks, but a user cannot carry SPARK's `tau`, kernel weights, or combined p-value decomposition across.

The speed comparison above is against SPARK-G (the Gaussian branch) and SpatialDE, not against SPARK's Poisson branch, because SPARK itself was excluded above 3,000 locations and did not converge on sparse simulations. The paper's own claim against SPARK is "orders of magnitude computational savings" in the abstract, with the concrete per-method numbers reported for the Gaussian/SpatialDE variants.

The test assumes a linear/projection kernel in coordinates plus Gaussian/cosine transformed coordinates; it is not a drop-in replacement for kernels outside that family, and the default `option="mixture"` runs 11 kernels (11x the statistic and p-value work of `option="single"`), so the cheapest configuration changes the output set.

SPARK-X is not exposed by omicverse 2.4.3; any omicverse comparison must call the R package or a reimplementation directly. Its `omicverse.space.svg` alternatives (Moran, SOMDE, SpatialDE) are different statistics and cannot be used to reproduce these numbers.

## Reuse plan for Bio-Opt

- ** analytic reformulation first**: ask whether the fitted model's test statistic has a closed-form null (weighted chi-square, Davies/Liu) that avoids estimating nuisance parameters. If yes, prove the statistic identity on a small matrix before timing.
- ** reuse shared state**: factor the `d x d` coordinate/kernel object once and reuse it across all genes and kernels; this is the step that makes cost linear in genes rather than per-gene iterative.
- ** prune unused work**: never materialise the `n x n` kernel; drop all-zero cells and genes once.
- ** vectorise/compile**: express the per-gene statistic as matrix products and parallelise the p-value map; keep the Davies call out of the inner data loop.
- **Verification**: compare detected gene sets and p-value ranks, not p-values, against SPARK on a small dataset where both converge; run a permutation null to check calibration before accepting the speed.

Pattern codes as used here (the referenced `pattern-catalog.md` was absent from this checkout):  prune unused work,  reuse shared/target-invariant state,  vectorise/compile,  analytic reformulation.
