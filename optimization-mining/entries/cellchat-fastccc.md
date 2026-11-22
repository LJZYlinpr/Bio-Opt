---
slug: cellchat-fastccc
method: CellChat
omicverse_binding: "omicverse.single.run_cellchat (native Python CellChat reimplementation) and omicverse.single.run_fastccc / ov.single.run_ccc(method='fastccc'); both exported in omicverse v2.4.3"
original:
  paper: "Jin et al., Nat Commun 12, 1088 (2021), doi:10.1038/s41467-021-21246-9"
  repo: "https://github.com/sqjin/CellChat"
  commit: "e4f68625b074247d619c2e488d33970cc531e17c"
optimized:
  paper: "Hou, Ma & Zhou, Nat Commun 16, 11428 (2025), doi:10.1038/s41467-025-66272-z"
  repo: "https://github.com/Svvord/FastCCC"
  commit: "cb6b9f3fc4c9b61a8c52fdcb2335ce2139926650"
optimization_type: algorithmic
mechanism: "Replace the label-permutation null with an analytic null for the communication score, obtained by FFT convolution of discretized single-cell expression distributions."
reported_speedup: "No CellChat-specific ratio is reported. FastCCC is at least 10x faster than nearly all other methods across 8 smaller datasets and is the only method that completes 1-2M-cell atlases, while CellChat (default parameters) hit memory errors at 100-200k cells (paper Fig. 3; HPC nodes with Intel Xeon E5-2683 v3, 500 GB RAM per node). Against CellPhoneDB, >=40x faster above 100k cells and ~500x on 2M cells, <20 min (paper Fig. 2g; Ubuntu 20.04, Intel Xeon E5-2620 v2, 250 GB RAM)."
equivalence: changed-output
optimization_families: [O1, O2, O6, O9]
sources:
  - "https://doi.org/10.1038/s41467-021-21246-9"
  - "https://github.com/sqjin/CellChat/blob/master/R/modeling.R"
  - "https://github.com/Svvord/FastCCC"
  - "https://github.com/Svvord/FastCCC/blob/main/README.md"
  - "https://github.com/Svvord/FastCCC/blob/main/fastccc/distrib.py"
  - "https://doi.org/10.1038/s41467-025-66272-z"
verified: true
---

## Baseline bottleneck

CellChat scores a ligand–receptor pair between two groups with a mass-action/Hill expression of the group-level ligand and receptor averages (CellChat paper, Methods; `computeCommunProb`, `R/modeling.R` L50–52, default `type = "triMean"`, `nboot = 100`, `Kh = 0.5`). Significance is a one-sided permutation test over cell-group labels: `M = 100` permutations by default (paper Eq. 3, https://doi.org/10.1038/s41467-021-21246-9), so the smallest attainable p-value is 0.01. The pinned commit is the current 2.x master, which retains the same permutation test (`nboot = 100`).

The cost is multiplicative in permutations and interactions. One permuted group-mean matrix is built up front, but for every ligand–receptor pair the code re-evaluates the per-group ligand/receptor summaries and the Hill score for all `nboot` permutations (`R/modeling.R` L226–250), then counts exceedances:

```r
permutation <- replicate(nboot, sample.int(nC, size = nC))
data.use.avg.boot <- my.sapply(X = 1:nboot, FUN = function(nE) { ... })
Pboot <- sapply(X = 1:nboot, FUN = function(nE) {
  data.use.avgB <- data.use.avg.boot[[nE]]
  dataLavgB <- computeExpr_LR(geneL[i], data.use.avgB, complex_input)
  dataRavgB <- computeExpr_LR(geneR[i], data.use.avgB, complex_input)
  dataLRB = Matrix::crossprod(dataLavgB, dataRavgB)
  .boot <- dataLRB^n/(Kh^n + dataLRB^n)
  ...
  return(as.vector(Pboot))
})
nReject <- rowSums(Pboot - Pnull > 0)
p = nReject/nboot
```

The repeated scoring is `O(nboot * nLR * nGroups^2)`, and the p-value grid is quantized by `nboot`. The paper reports that CellChat runs out of memory on 100–200k cells with default parameters (FastCCC paper, Fig. 3 discussion), so the permutation design also blocks atlas-scale use.

## Optimization trajectory

1. **Keep score families, replace the null by an analytic distribution .** FastCCC derives the null law of the score under random group labels instead of sampling labels. Each gene's expression inside a group is discretized into a pmf (0.01-wide bins after log normalization); the distribution of a sum over `n` cells is the `n`-fold convolution of that pmf, computed with FFT (`fastccc/distrib.py`, `Distribution.__add__` / `__pow__`). Cost per interaction is a small number of convolutions, with no `nboot` factor and no 1/100 p-value floor (paper, "A general CCC testing framework").
2. **Exact order-statistic nulls for quantile scores .** The median/Q3/Q90 variants score an order statistic rather than a mean. FastCCC computes the exact pmf of the k-th order statistic from the same discrete distribution, including ties, instead of resampling labels (paper Eq. 5 and following; `fastccc/dist_complex.py`, `dist_iid_set.py`).
3. **Normal branch for common cell types .** When the group has `>= 30` cells the sum is approximated by a normal with the exact mean and variance, skipping the FFT convolution (`distrib.py`, `__pow__`: `if n >= 30 ... return Distribution(dtype='normal', ...)`). The approximation is used only where the central limit applies; rare groups keep the exact discrete law.
4. **Batch the 16 score variants and combine p-values .** Four expression summaries (mean, median, Q3, Q90) x two complex aggregations (minimum, average) x two L-R combinations (arithmetic, geometric) are computed from the same per-gene pmfs, and the 16 p-values are combined with the Cauchy combination test (paper, "Cauchy combination"; `fastccc/cauchy_combine.py`). More statistics per pass, still no permutations.
5. **Reuse reference nulls (reuse target-invariant state).** For reference-based analysis, per-gene distributions of a healthy panel are stored and reused, so the query is tested against precomputed nulls rather than re-deriving them (README, "Reference panels"; `fastccc/build_reference.py`).

```python
def __add__(self, pmf2):
    ...
    new_pmf_array = fftconvolve(self.get_pmf_array, pmf2.get_pmf_array)
    ...
def __pow__(self, n):
    if n >= 30 or self.is_normal_type:
        return Distribution(dtype='normal', loc=self.get_mean * n,
                            scale=np.sqrt(n) * self.get_std)
    new_pmf_array = self.get_pmf_array
    for i in range(1, n):
        new_pmf_array = fftconvolve(new_pmf_array, self.get_pmf_array)
```

Why this removes work: the permutation loop disappears, so runtime no longer scales with `nboot` and p-values are not discretized to `1/nboot`. What remains is per-gene discretization plus convolutions. The paper's own accounting shows the computation is no longer the bottleneck: after this change about 25% of runtime is spent reading data and about 60% of memory is the input matrix itself (FastCCC paper, Fig. 2g,h discussion).

## Evidence

- FastCCC (https://doi.org/10.1038/s41467-025-66272-z): "FastCCC is at least 40 times faster than CPDB when the number of cells exceeds 100,000 and nearly 500 times faster on a dataset with two million cells, completing all calculations under 20 minutes in the latter case" (paper, Results). Hardware for that comparison: Ubuntu 20.04, Intel Xeon E5-2620 v2, 250 GB RAM (paper, Methods).
- Against 7 methods on 8 datasets: "FastCCC was at least ten times faster than nearly all other methods, while also consuming the least memory"; "CellChat and CellCall encountered memory errors when processing datasets with 100–200 k cells using default parameters" (paper, Fig. 3 discussion, https://doi.org/10.1038/s41467-025-66272-z; Slurm cluster, Intel Xeon E5-2683 v3, 500 GB per node).
- Accuracy: "CPDB requires at least 10,000 permutations—ten times its default setting—to match FastCCC's p-value accuracy" (paper, Results, Fig. 2f).
- Runtime/memory accounting: ~60% of memory and ~25% of compute time are data reading (paper, Fig. 2g,h).

## Equivalence and caveats

`changed-output`. FastCCC is not a port of CellChat's statistic: CellChat's probability is a Hill-normalized product with agonist, antagonist and co-receptor modulation and `triMean` group averages. FastCCC's 16 modules are built from mean/median/Q3/Q90 summaries, minimum/average complex aggregation and arithmetic/geometric L-R combination; none includes the Hill or cofactor terms, and its single mean/minimum/arithmetic choice reproduces the CellPhoneDB score, not CellChat's (README scoring table; paper, "A general CCC testing framework"). Communication probabilities, the tested interaction set, and therefore the significant network differ. FastCCC's own equivalence evidence is against CellPhoneDB, not CellChat: p-values are "nearly identical" to CPDB once CPDB uses 10,000 permutations (paper, Fig. 2f).

Further caveats: FastCCC's p-values assume independence between ligand and receptor expression summaries under the null (paper, Methods); the analytic tail below the discretization floor is approximated; output is tables + HTML, not CellChat's R network objects, so CellChat's downstream centrality/pattern functions are not interchangeable. Re-validate on a matched LR database and matched summary statistic before comparing numbers.

## omicverse binding

`omicverse.single.run_cellchat` is a native Python reimplementation of CellChat 2.x defaults (`n_perms=100`, Hill/agonist/antagonist/co-receptor terms, `omicverse/external/cellchat/_engine.py`), so it inherits the baseline bottleneck. `omicverse.single.run_fastccc` calls a vendored engine (`omicverse/external/fastccc`) with CellPhoneDB v5 interactions and ignores `n_perms` because p-values are analytic; `ov.single.run_ccc(method='fastccc')` dispatches to the same backend. Both symbols are exported in omicverse v2.4.3 (`omicverse/single/__init__.py` L65–84, L370–391).

## Reuse plan for Bio-Opt

Primary pattern: ** analytic reformulation** — derive the null law of the test statistic instead of resampling. Secondary: ** vectorize/compile** the convolution and ** controlled approximation** (normal branch, fixed 0.01 discretization) with an explicit error check.

Test first: on one small dataset, compute CellChat's observed probability and compare its 100-permutation p-value with an FFT-convolution null built for the *same* product score (no Hill), sweeping the discretization width. Accept if the analytic tail matches a 10,000-permutation reference within Monte-Carlo error across the p-value range; then measure runtime versus `nboot`. Do not transfer FastCCC p-values to CellChat's Hill score without re-deriving the null.
