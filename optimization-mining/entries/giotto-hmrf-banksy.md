---
slug: giotto-hmrf-banksy
method: "Giotto HMRF -> BANKSY (spatial domain segmentation)"
omicverse_binding: "omicverse.space.clusters(methods=['Banksy']) imports the Python pybanksy package (Banksy_py) and calls initialize_banksy + run_banksy_multiparam (omicverse/space/_cluster.py:651-692); Giotto HMRF has no omicverse binding in the inspected 2.4.3 tree (repository-wide grep found no Giotto/HMRF reference)"
original:
  paper: "Dries R et al. Giotto: a toolbox for integrative analysis and visualization of spatial expression data. Genome Biol 22:78 (2021), doi:10.1186/s13059-021-02286-2"
  repo: "https://github.com/drieslab/Giotto"
  commit: "b4a774fc9d7a83a326991949d2ef2f637c139872"
optimized:
  paper: "Singhal V et al. BANKSY unifies cell typing and tissue domain segmentation for scalable spatial omics data analysis. Nat Genet 56:431-441 (2024), doi:10.1038/s41588-024-01664-3"
  repo: "https://github.com/prabhakarlab/BANKSY"
  commit: "70175fbfe10d5ee38efeaba86575bc55965c1ed3"
optimization_type: algorithmic
mechanism: "Replace the Potts-model HMRF label field (50 beta values, each refitting the label field in an external smfishHmrf process) with a single feature augmentation (own expression concatenated with kNN neighborhood means), one PCA, and graph clustering on the resulting embedding."
reported_speedup: "10-1,000x faster and 2-60x more scalable than existing spatial clustering methods, including Giotto's HMRF module; BANKSY was the only method that scaled to 2 million MERFISH cells, all methods benchmarked on a 16-thread, 128 GB machine with a 16 h cutoff (Nat Genet 2024, Fig. 6 and Results)"
equivalence: changed-output
optimization_families: [O1, O5, O6, O11]
sources:
  - "https://doi.org/10.1186/s13059-021-02286-2"
  - "https://github.com/drieslab/Giotto/blob/b4a774fc9d7a83a326991949d2ef2f637c139872/R/python_hmrf.R"
  - "https://doi.org/10.1038/s41588-024-01664-3"
  - "https://pmc.ncbi.nlm.nih.gov/articles/PMC10937399/"
  - "https://github.com/prabhakarlab/BANKSY"
  - "https://github.com/prabhakarlab/Banksy_py/blob/9278996c39e376277d57ef95278000447ba6c57c/src/banksy/main.py"
  - "https://github.com/prabhakarlab/Banksy_py/blob/9278996c39e376277d57ef95278000447ba6c57c/src/banksy/embed_banksy.py"
  - "https://github.com/Starlitnightly/omicverse/blob/d0c3325ca72404e28ec5161a2f71ea82db531974/omicverse/space/_cluster.py"
  - "https://www.biorxiv.org/content/10.1101/2023.01.10.523386v2"
verified: true
---

## Baseline bottleneck

This entry replaces the requested BayesSpace -> CellCharter pair. CellCharter's own benchmark on 12 joint DLPFC samples reports that "its computational efficiency was second only to BayesSpace, which was slightly faster although requiring more memory" (CellCharter preprint, Results, Fig. 1c; Nat Genet 2024 version keeps the same comparison), so CellCharter is a memory/scalability optimisation, not a speed trajectory. BANKSY instead reports a measured 10-1,000x speedup against the HMRF family, including Giotto's HMRF module, which is why the HMRF -> BANKSY pair is used here.

Giotto's `doHMRF` fits a hidden Markov random field with a Potts prior on the spatial label field. Two structural costs dominate:

- **Beta sweep inside one external process.** The default `betas = c(0, 2, 50)` means start 0, increment 2, 50 values. Giotto writes cell locations, spatial genes, the spatial network and the expression matrix to CSV, then calls `system` to launch the `smfishHmrf` Python reader with all beta values (`R/python_hmrf.R:54`, `R/python_hmrf.R:296-320`). Fresh label-field iterations are run for every beta, and the best beta is selected afterwards: 50 full HMRF fits per dataset plus subprocess and disk round-trips.
- **Iterative label sampling.** Within each run, HMRF updates labels cell by cell over neighbouring cells for many sweeps; BayesSpace, the other HMRF implementation in the same benchmark, documents 50,000 MCMC iterations with a 1,000-iteration burn-in as its default for Visium (BANKSY Methods, citing the BayesSpace vignette).

The work is therefore proportional to `n_beta x n_sweeps x n_cells`, with the beta dimension and the label field both removed by the alternative below.

## Optimization trajectory

1. **Replace the label field with augmented features .** BANKSY builds a kNN graph in space, computes each cell's neighbourhood mean expression, z-scores each block, weights the neighbourhood blocks by `lambda`, and concatenates them with the cell's own expression. The output is one `n x (G * (m+1))` matrix; there is no latent label field, no Potts prior, no beta, and no sampling (`Banksy_py/src/banksy/embed_banksy.py:25-97`, `Banksy_py/src/banksy/main.py:313-350`).
2. **One embedding, one clustering pass .** PCA is run on the augmented matrix (20 PCs in the Python implementation) and Leiden/Louvain clusters that embedding, so the whole domain segmentation is one deterministic pipeline instead of `n_beta` HMRF fits (`cluster_methods.py`; omicverse calls `run_banksy_multiparam`).
3. **Share the neighbourhood work across lambdas .** `create_nbr_matrix` and `create_mat_list` are computed once per decay scheme; the multi-parameter runner then forms each lambda's augmented matrix from the same cached list (`embed_banksy.py:41-66`, `main.py:313-350`). A lambda sweep therefore repeats only z-scoring/concatenation/PCA/clustering, not graph construction.
4. **Bound memory with a lazy operator .** The R implementation's `runBanksyPCA` computes PCA without materialising the full augmented matrix: the kNN graph and BANKSY transform are applied on the fly inside the iterative PCA solver, which the README describes as "memory-efficient and scales to millions of cells". The Python implementation used by omicverse instead densifies (`mat.todense`), which is the main caveat for the pybanksy path.
5. **Single hyperparameter switches task.** `lambda` near 0.2 targets cell typing, near 0.8 domain segmentation; no model refit per beta is needed.

```python
# Banksy_py/src/banksy/main.py:341-350 (concatenate_all)
scaled_list = []
for n in range(len(matrix_list)):
    mat = matrix_list[n]
    if issparse(mat):
        mat = mat.todense
    if np.iscomplexobj(mat):
        mat = np.absolute(mat)
    scaled_list.append(scale_factors[n] * zscore(mat, axis=0))
concatenated_matrix = np.concatenate(scaled_list, axis=1)
```

The work it replaces (Giotto R/python_hmrf.R:296-320):

```r
betas_param <- c("-b", betas)
betas_final <- paste(betas_param, collapse = " ")
reader_command <- paste0(python_path, " ", reader_path,
    " -l ", cell_location, " -g ", spatial_genes,
    " -n ", spatial_network, " -e ", expression_data,
    " -o ", output_data, " -a ", name, " -k ", k,
    " ", betas_final, " -t ", tolerance,
    " -z ", zscore, " -s ", seed, " -i ", numinit)
system(command = reader_command)
```

## Evidence

- "In summary, BANKSY is 10-1,000 times faster and 2-60 times more scalable than existing spatial clustering algorithms." The runtime benchmark includes "Giotto's HMRF module" and BayesSpace, up to 2 million cells, on a "virtual machine with 128 GB and 16 threads (Ubuntu v.18.04.5 LTS)", with runs terminated at 16 h (Nat Genet 2024, Results, Fig. 6 and Methods, Scalability analysis).
- "BANKSY was the only method that scaled to 2 million cells, while the other spatial algorithms hit memory or software limits, or failed to return a result within 16 h"; "BANKSY's runtime was only marginally higher than that of nonspatial clustering, while the other methods were 1-3 orders of magnitude slower."
- Accuracy, 12 DLPFC Visium sections: "By all three metrics, BANKSY outperformed Giotto, BayesSpace, SpaGCN and SpiceMix" (ARI, NMI, MCC; Fig. 5a,b, Extended Data Fig. 6). The comparison used Giotto v1.1.0 with beta swept over the authors' recommended grid in steps of 2 (the best/stable beta selected per dataset) and BayesSpace v1.5.1 with its 50,000-iteration default, i.e. a deliberately strong HMRF baseline.
- CellCharter counter-evidence: on 12 joint DLPFC samples "its computational efficiency was second only to BayesSpace, which was slightly faster although requiring more memory" (biorxiv v2, Results) - hence the pair substitution documented above.

## Equivalence and caveats

`changed-output`: BANKSY does not fit an HMRF and produces no posterior probabilities, no beta, and no per-domain variance parameters. Cluster labels are expected to differ from Giotto HMRF even when the number of clusters matches; equivalence must be judged against annotation (ARI/NMI/MCC) or marker enrichment, not label identity.

Two BANKSY parameters change the scientific target: `lambda` (0 for non-spatial, ~0.2 for cell typing, ~0.8 for domain segmentation) and `k_geom` (spatial neighbours). Giotto's HMRF has no comparable cell-typing/domain-segmentation knob, so a fair comparison must fix the number of domains and tune both methods on the same reference.

The 10-1,000x claim was measured on the implementation used in the paper's benchmark (the R/Bioconductor pipeline, timed with peakRAM). The Python package used by omicverse densifies the augmented matrix (`mat.todense` in `concatenate_all`), so the million-cell scaling and peak memory should not be assumed for `ov.space.clusters(methods=['Banksy'])` without measurement.

Giotto also provides HMRF-based visualization and network-based downstream steps that BANKSY does not replace; only the domain-segmentation step is compared here.

## Reuse plan for Bio-Opt

omicverse binding: `ov.space.clusters(adata, methods=['Banksy'], methods_kwargs={'Banksy': {...}})` (`omicverse/space/_cluster.py:651`) requires `pip install pybanksy`; it stores per-parameter embeddings as `adata.obsm['X_banksy_*']`. There is no omicverse wrapper for Giotto HMRF, so the baseline must be run from R or a separate environment.

- ** analytic reformulation first**: when the expensive object is a sampled label field, check whether the same spatial information can be encoded as deterministic features (neighbourhood moments) plus a standard clustering step; validate with ARI against annotation.
- ** prune unused work**: remove the beta sweep/selection loop; it exists only because the HMRF has a regularization parameter, which the feature-augmentation route does not.
- ** reuse shared state**: build the kNN graph and neighbourhood matrices once, then reuse them across lambda/resolution sweeps.
- ** bound memory**: use the lazy PCA operator (or chunked/sparse augmentation) before scaling to millions of cells; do not port the dense Python path to large data.
- **Verification**: compare domain labels by ARI/NMI/MCC against the same annotation and fix the number of domains; also check that the chosen lambda is reported, since it changes the cell-typing vs domain-segmentation regime.

Pattern codes as used here (the referenced `pattern-catalog.md` was absent from this checkout):  prune unused work,  reuse shared/target-invariant state,  stream/chunk to bound memory,  analytic reformulation.
