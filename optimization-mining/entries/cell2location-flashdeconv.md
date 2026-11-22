---
slug: cell2location-flashdeconv
method: cell2location
omicverse_binding: "omicverse.space.Deconvolution(method='cell2location'|'FlashDeconv'); both branches present in omicverse v2.4.3"
original:
  paper: "Kleshchevnikov et al., Nat Biotechnol 40, 661-671 (2022), doi:10.1038/s41587-021-01139-4"
  repo: "https://github.com/BayraktarLab/cell2location"
  commit: "afb07fda7f89458e44a4afdf21c26aeb7251ebfa"
optimized:
  paper: "Yang, Chen & Zhang, bioRxiv 2025, doi:10.64898/2025.12.22.696108"
  repo: "https://github.com/cafferychen777/flashdeconv"
  commit: "a7c90177beb1c5aacf368e91e81167ab4381349a"
optimization_type: approximation
mechanism: "Replace hierarchical Bayesian negative-binomial inference with a deterministic leverage-weighted CountSketch plus a graph-regularized non-negative least-squares fit solved by a parallel block coordinate descent."
reported_speedup: "Not reported as a direct cell2location ratio. FlashDeconv: 1,000,000 Visium HD bins in 40 s and 12.6 GB peak on 32 CPU threads; 100,000 bins 5.4 s / 2.7 GB; 10,000 bins 1.9 s / 1.9 GB (README, colorectal-cancer Visium HD, 18,082 genes x 38 cell types). The preprint abstract claims 'accelerating inference by orders of magnitude' versus Bayesian deconvolution; cell2location's own benchmark reports 8.66 min on GPU for its simulated case (preprint, Results)."
equivalence: approximate
optimization_families: [O1, O3, O5, O6, O9, O10, O11]
sources:
  - "https://doi.org/10.1038/s41587-021-01139-4"
  - "https://doi.org/10.1101/2020.11.15.378125"
  - "https://github.com/BayraktarLab/cell2location"
  - "https://github.com/BayraktarLab/cell2location/blob/master/docs/notebooks/cell2location_tutorial.ipynb"
  - "https://github.com/BayraktarLab/cell2location/blob/master/cell2location/models/_cell2location_module.py"
  - "https://github.com/cafferychen777/flashdeconv"
  - "https://github.com/cafferychen777/flashdeconv/blob/main/README.md"
  - "https://github.com/cafferychen777/flashdeconv/blob/main/flashdeconv/core/solver.py"
  - "https://github.com/cafferychen777/flashdeconv/blob/main/flashdeconv/core/sketching.py"
verified: partial
---

## Baseline bottleneck

cell2location is a hierarchical Bayesian model. A reference model learns cell-type signatures, then a spatial model estimates a non-negative abundance `w_sf` for every cell type at every location, with a Gamma prior whose mean is a factorization across cell types (to borrow strength between similar locations) and a Gamma–Poisson likelihood that absorbs technical effects (`cell2location/models/_cell2location_module.py` L354–372, L461–476, https://github.com/BayraktarLab/cell2location/blob/master/cell2location/models/_cell2location_module.py):

```python
with obs_plate as ind:
    w_sf_mu = z_sr_groups_factors @ x_fr_group2fact
    w_sf = pyro.sample("w_sf", dist.Gamma(w_sf_mu * self.w_sf_mean_var_ratio_tensor,
        self.w_sf_mean_var_ratio_tensor))   # (n_obs, n_factors)
...
mu = ((w_sf @ self.cell_state) * m_g + (obs2sample @ s_g_gene_add)) * detection_y_s
with obs_plate:
    pyro.sample("data_target",
                dist.GammaPoisson(concentration=alpha, rate=alpha / mu),
                obs=x_data)
```

The observation plate has one latent vector per location (`create_plates`, `pyro.plate("obs_plate", size=n_obs)`), so variational training touches every location in every epoch. The official tutorial trains the reference model for 250 epochs and the spatial model for 30,000 epochs with `train_size=1`, then draws 1,000 posterior samples, all on GPU. Cost therefore scales with locations x epochs x posterior samples, and the deliverable is a posterior (q05/q50/q95), not a point estimate. The README lists "scalability to 100k-mln+ locations using amortised inference" as a future direction, which indicates that per-location inference is the current scaling limit.

## Optimization trajectory

1. **Replace the generative model by a penalized regression .** FlashDeconv estimates non-negative abundances `beta` by minimizing `0.5*||Y_w - beta*X_w||^2 + 0.5*lambda*Tr(beta^T L beta) + rho*||beta||_1`, then row-normalizes to proportions (README, "How it works"). This drops the likelihood, the variational training loop, and posterior sampling; the fit is a deterministic optimization with a fixed iteration cap.
2. **Reduce genes before fitting .** It keeps highly variable genes of the spatial data plus marker genes of each reference type, weights each gene by its leverage score in the reference (how strongly it separates cell types), and log-normalizes both matrices (README steps 1–2; `flashdeconv/core/preprocessing.py`, `utils/genes.py`). Discriminative genes dominate and uninformative genes are not carried through every iteration.
3. **Sketch the gene axis .** A leverage-weighted CountSketch hashes genes into `sketch_dim` buckets with random signs and amplitudes scaled by `sqrt(leverage)` (`flashdeconv/core/sketching.py` `build_countsketch_matrix`, default 512). Each regression then works in 512 dimensions instead of ~18k genes, while preserving rare-type signal in expectation. The sketch is deterministic for a fixed seed.
4. **Precompute the normal equations .** `bcd_solve` computes the Gram matrix `X_sketch @ X_sketch.T` (K x K, K = number of cell types) and `H = X_sketch @ Y_sketch.T` once (`solver.py` L360–370). Per-iteration cost becomes O(n_bins * K^2) with K ~ 38, rather than O(n_bins * n_genes * K).
5. **Solve with a fused, parallel block coordinate descent .** `_bcd_iteration_fused` loops over bins with Numba `prange`, updates all cell types for a bin from its precomputed `H[:, i]`, adds the kNN-graph neighbor mean (Laplacian penalty), applies soft-thresholding (L1) and non-negativity, and computes convergence statistics in the same pass with double buffering (`solver.py` L113–196). The solver is linear in bins and runs in Numba's native parallel mode without Python object overhead.
6. **Replace the hierarchical prior by a kNN graph penalty .** Spatial sharing between neighboring bins is a graph Laplacian term instead of a factorized Gamma prior, so there is no per-location posterior to sample; `lambda_spatial` is scaled automatically (README).

```python
@jit(nopython=True, parallel=True, cache=True)
def _bcd_iteration_fused(H, XtX, beta_in, beta_out, neighbor_indices,
                         neighbor_indptr, lambda_, rho, spot_diffs, spot_abs):
    for i in prange(n_spots):
        Xty_i = H[:, i]
        ...
        neighbor_sum = np.zeros(n_cell_types)
        for idx in range(start, end):
            j = neighbor_indices[idx]
            for k in range(n_cell_types):
                neighbor_sum[k] += beta_in[j, k]
        update_spot_with_Xty(Xty_i, XtX, beta_out[i], neighbor_sum,
                             n_neighbors, lambda_, rho)
```

Why this removes work: the 30,000-epoch GPU training and 1,000-sample posterior are replaced by a fixed number of BCD sweeps over a 512-dimensional sketch, with the expensive matrix products precomputed once. Accuracy is retained by weighting and selecting genes rather than by a full generative likelihood; the README reports mean Pearson 0.946 on Spotless versus 0.918 for cell2location.

## Evidence

- README performance table (https://github.com/cafferychen777/flashdeconv#performance; Visium HD colorectal cancer, 18,082 genes, 38 cell types, 32 CPU threads): 10k bins 1.9 s / 1.9 GB; 100k bins 5.4 s / 2.7 GB; 1M bins 40 s / 12.6 GB.
- README accuracy on the 54 silver-standard Spotless datasets (mean Pearson with true proportions): FlashDeconv 0.946, RCTD 0.934, cell2location 0.918.
- Preprint abstract: "Benchmarking demonstrates accuracy comparable to top-tier Bayesian methods while accelerating inference by orders of magnitude"; applied to ovarian cancer cohorts it "reproduces clinical response signatures in seconds".
- Baseline cost: cell2location tutorial trains 250 epochs (reference) + 30,000 epochs (spatial) on GPU and draws 1,000 posterior samples; its preprint reports 8.66 min on GPU versus Stereoscope's 14.8 h on the same simulated case.

## Equivalence and caveats

`approximate`. Both methods estimate cell-type composition from a single-cell reference, and FlashDeconv's proportions are close to cell2location's on the Spotless benchmark, but the estimands differ: cell2location returns absolute cell abundance with a posterior (q05/q50/q95) and can export cell-type-specific imputed expression, while FlashDeconv returns relative proportions summing to 1 plus a dominant type, with model-based (not posterior) uncertainty. The spatial prior also differs (graph Laplacian vs hierarchical Gamma), so per-spot values are not expected to be identical.

Practical caveats: the reference must contain every cell type present, otherwise the fit misassigns transcripts — FlashDeconv ships a reference diagnostic that flags unexplainable bins and names missing types (README); the sketch is a randomized approximation, so check sensitivity to `sketch_dim`; `lambda_spatial` default is `"auto"` in the package, but omicverse passes `lambda_spatial=5000.0` and `sketch_dim=512` explicitly, which is a different operating point from the README's defaults. For segmented platforms (Xenium) the README recommends `lambda_spatial=0`.

## omicverse binding

`omicverse.space.Deconvolution` exposes `method='cell2location'` (vendored `omicverse/external/space/cell2location`, default `max_epochs=250`/`30000`, GPU defaults applied) and `method='FlashDeconv'` (imports the `flashdeconv` package and calls `fd.tl.deconvolve`, storing `obsm['flashdeconv']`). Both branches and the `flashdeconv_kwargs` parameter are present in v2.4.3 (`omicverse/space/_deconvolution.py` L252–305, L350–355, L401, L454, L483–512), so the two estimators can be swapped behind one API.

## Reuse plan for Bio-Opt

Patterns: ** controlled approximation** (leverage-weighted CountSketch), ** compile/fuse**, ** parallelize**, ** prune genes**.

Test first: on one Visium slide, compare FlashDeconv proportions with cell2location q50 (Pearson per cell type, and rank of dominant type), then sweep `sketch_dim` in {128, 256, 512, 1024} and `lambda_spatial` to find where accuracy saturates. Accept the sketch if accuracy loss is within the pre-declared tolerance while runtime scales with `sketch_dim`; keep cell2location as the reference for absolute abundance and posterior intervals, which FlashDeconv does not provide.
