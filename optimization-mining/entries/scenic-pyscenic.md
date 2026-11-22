---
slug: scenic-pyscenic
method: "SCENIC (R, GENIE3) -> pySCENIC (GRNBoost2 + distributed per-target fits)"
omicverse_binding: "omicverse.single.SCENIC and omicverse.single.grn (vendored arboreto: grnboost2/genie3, and vendored ctxcore FeatherRankingDatabase)"
original:
  paper: "Aibar et al., Nature Methods 2017, doi:10.1038/nmeth.4463"
  repo: "https://github.com/aertslab/SCENIC"
  commit: "7a74341745cecd3505310c6c5755cad456756cf9"
optimized:
  paper: "Van de Sande et al., Nature Protocols 2020, doi:10.1038/s41596-020-0336-2; GRNBoost2: Moerman et al., Bioinformatics 2019, doi:10.1093/bioinformatics/bty916"
  repo: "https://github.com/aertslab/pySCENIC"
  commit: "pySCENIC 06bafba412792f6efa5a552a23bb221cc3bdea1b; arboreto 79f916b0ea25c00989331b8db243826049c3d66c"
optimization_type: "algorithmic + system"
mechanism: "Replace per-target random-forest GENIE3 with early-stopped stochastic gradient boosting (GRNBoost2) over a fixed TF predictor matrix, and run the independent per-target regressions through a process pool or dask instead of a serial R loop."
reported_speedup: "10x versus R SCENIC (Nat Protoc 2020 abstract: pySCENIC 'resulting in a tenfold increase in speed'); '<2 h' for 10,000 genes x 50,000 cells (same abstract). Per-step GRNBoost2-versus-GENIE3 factor: not reported in the accessible abstract"
equivalence: changed-output
optimization_families: [O8, O9, O10, O11]
sources:
  - "https://doi.org/10.1038/s41596-020-0336-2"
  - "https://doi.org/10.1093/bioinformatics/bty916"
  - "https://doi.org/10.1038/s41592-023-01938-4"
  - "https://github.com/aertslab/SCENIC/blob/master/README.md"
  - "https://github.com/aertslab/SCENIC/blob/master/R/runGenie3.R"
  - "https://github.com/aertslab/pySCENIC/blob/06bafba412792f6efa5a552a23bb221cc3bdea1b/README.rst"
  - "https://github.com/aertslab/pySCENIC/blob/06bafba412792f6efa5a552a23bb221cc3bdea1b/src/pyscenic/cli/arboreto_with_multiprocessing.py"
  - "https://github.com/aertslab/arboreto/blob/79f916b0ea25c00989331b8db243826049c3d66c/arboreto/core.py"
  - "https://github.com/aertslab/arboreto/blob/79f916b0ea25c00989331b8db243826049c3d66c/README.rst"
  - "https://github.com/aertslab/SCENICplus/blob/main/Dockerfile"
verified: true
---

## Scope note: SCENIC+ is not the speed-optimized successor

The requested pair `pySCENIC -> SCENIC+` does not hold as a performance trajectory. SCENIC+ (Bravo Gonzalez-Blas et al., Nat Methods 2023) is an *enhancer-driven* GRN method for joint scRNA-seq + scATAC-seq: it predicts genomic enhancers, links them to target genes and clusters >30,000 motifs. Its abstract and README make no runtime claim against pySCENIC, and its pipeline adds pycisTopic topic modelling with MALLET (installed explicitly in its Dockerfile) plus pycistarget, i.e. strictly more work on more modalities — https://doi.org/10.1038/s41592-023-01938-4, https://github.com/aertslab/scenicplus, and https://github.com/aertslab/SCENICplus/blob/main/Dockerfile (pycisTopic + /opt/mallet). The SCENIC+ repository does carry `BASCA_speedup` and `BASCA_speedup2` branches, which confirms that speed work is a separate, non-mainline effort rather than a property of the released SCENIC+.

The actual accelerated successor of SCENIC is pySCENIC, and the R repository says so directly: "SCENIC is deprecated, use pySCENIC instead" (https://github.com/aertslab/SCENIC/blob/master/README.md). The pySCENIC README names the mechanism: "The most impactful speed improvement is introduced by the arboreto package in step 1. This package provides an alternative to GENIE3 called GRNBoost2" (https://github.com/aertslab/pySCENIC/blob/06bafba412792f6efa5a552a23bb221cc3bdea1b/README.rst). This entry therefore records `SCENIC (R) -> pySCENIC`, which is a genuine execution trajectory.

## Baseline bottleneck

`runGenie3` in R SCENIC calls `GENIE3::GENIE3(exprMat, regulators = inputTFs, nCores = nCores, targets = genesSplit[[i]], ...)` one target-part at a time (`R/runGenie3.R:90`). GENIE3 fits a random-forest regression per target gene, so work is `#targets x #trees x (#samples x #TFs)` split searches, executed mostly serially and materialized as a dense weight matrix per part. Arboreto's GENIE3 profile (the Python baseline that pySCENIC replaced) is explicit about that cost: `RF_KWARGS = {'n_jobs': 1, 'n_estimators': 1000, 'max_features': 'sqrt'}` (`arboreto/core.py:28-32`) — 1,000 trees per target.

The downstream steps add their own serial loops: regulon pruning against motif rankings (cisTarget) and AUCell scoring per regulon.

## Optimization trajectory

1. **Change the estimator, not just its parallelism** — GRNBoost2 keeps the GENIE3 architecture (one regression per target, TF-only predictors) but uses stochastic gradient boosting with feature subsampling and out-of-bag early stopping:

```python
# arboreto/core.py:49-55
SGBM_KWARGS = {
    'learning_rate': 0.01,
    'n_estimators': 5000,  # can be arbitrarily large
    'max_features': 0.1,
    'subsample': 0.9
}
# core.py:136
regressor.fit(tf_matrix, target_gene_expression,
              monitor=EarlyStopMonitor(early_stop_window_length))   # window = 25
```

Why less work: `max_features=0.1` cuts split search to 10% of candidate TFs per tree, `subsample=0.9` uses 90% of cells per tree, and early stopping means the 5,000-tree ceiling is never reached when the OOB improvement plateaus (controlled approximation). The defaults trade a small, measurable accuracy change for a much shorter training path.

2. **Build the predictor matrix once** — `to_tf_matrix` extracts the TF columns into one `n_cells x n_TFs` array that every target regression reads (reuse target-invariant state); only the target vector changes between tasks. This also shrinks each regression from "all genes as regulators" to "TFs only", which is the biologically intended model (prune unused work).

3. **Dispatch per-target fits as independent tasks** — each target's regression depends only on the shared TF matrix and its own expression vector, so the loop is embarrassingly parallel. pySCENIC ships a dedicated CLI that replaces dask with a `multiprocessing.Pool` and feeds target indices through `imap(chunksize=1)`:

```python
# pySCENIC/src/pyscenic/cli/arboreto_with_multiprocessing.py:168-186
with Pool(args.num_workers) as p:
    adjs = list(tqdm(p.imap(partial(run_infer_partial_network,
                gene_names=gene_names, ex_matrix=ex_matrix,
                tf_matrix=tf_matrix, tf_matrix_gene_names=tf_matrix_gene_names,
                method_params=method_params, seed=args.seed),
        target_gene_indices(gene_names, target_genes="all"),
        chunksize=1), total=len(gene_names)))
```

Why it can help: each per-target fit is small relative to dask's graph construction and scheduler round-trips, so a plain process pool removes that per-task overhead (release 0.12.0 made the custom pool the default; the note does not quantify the gain). dask remains available for multi-node runs (parallelism).

4. **Share the big matrix instead of shipping it ** — the AUCell step does what the GRN pool does not: it copies the ranking matrix once into a `multiprocessing.RawArray`, forks one child per signature chunk, and each child reads the same read-only block and writes into a preallocated shared result buffer:

```python
# pySCENIC/src/pyscenic/aucell.py:138-160
shared_ro_memory_array = RawArray(DTYPE_C, mul(*df_rnk.shape))
array = np.frombuffer(shared_ro_memory_array, dtype=DTYPE)
array[:] = df_rnk.values.ravel(order="C")          # one copy, not per task
auc_mtx = Array("d", len(cells) * len(signatures)) # shared result buffer
chunk_size = ceil(float(len(signatures)) / num_workers)
processes = [Process(target=_enrichment,
                     args=(shared_ro_memory_array, chunk, genes, cells,
                           auc_threshold, auc_mtx, (chunk_size * len(cells)) * idx))
             for idx, chunk in enumerate(chunked(signatures, chunk_size))]
```

   cisTarget pruning follows the same shape per TF/module over a precomputed feather ranking database (`memory_map=False`; faster feather loading added in 0.11).

5. **Cheaper I/O and memory** — rankings use `uint32`, sparse matrices are accepted by the GRN step (`--sparse`), and intermediate files are compressed; release 0.12.1 also replaced `flatten` with `ravel` in AUCell to avoid a copy (the `ravel(order="C")` above).

Relevant locations: `arboreto/core.py:28-60` (regressor profiles), `arboreto/core.py:86-160` (`to_tf_matrix`, `fit_model`, early stopping), `arboreto/core.py:281-330` (`infer_partial_network`), `pySCENIC/src/pyscenic/cli/arboreto_with_multiprocessing.py:101-186`, `pySCENIC/src/pyscenic/aucell.py:130-180` (shared-memory AUCell), `SCENIC/R/runGenie3.R:50-95` (per-part GENIE3 with checkpoint files).

## Evidence

- Nat Protoc 2020 abstract: "SCENIC has been refactored and reimplemented in Python (pySCENIC), resulting in a tenfold increase in speed... For a data set of 10,000 genes and 50,000 cells, the pipeline runs in <2 h." — https://doi.org/10.1038/s41596-020-0336-2
- GRNBoost2 abstract: "GRNBoost2 is an efficient algorithm for regulatory network inference using gradient boosting, based on the GENIE3 architecture. Arboreto is a computational framework that scales up GRN inference algorithms complying with this architecture." — https://doi.org/10.1093/bioinformatics/bty916
- Arboreto README: the library "addresses this issue by providing a computational strategy that allows executing the class of GRN inference algorithms exemplified by GENIE3 on hardware ranging from a single computer to a multi-node compute cluster"; GRNBoost2 is "a novel and fast GRN inference algorithm using Stochastic Gradient Boosting Machine regression with early-stopping regularization" — https://github.com/aertslab/arboreto/blob/79f916b0ea25c00989331b8db243826049c3d66c/README.rst
- pySCENIC README: dask "for distributed computing", and GRNBoost2 as "the most impactful speed improvement" — https://github.com/aertslab/pySCENIC/blob/06bafba412792f6efa5a552a23bb221cc3bdea1b/README.rst
- Hardware for the 10x and <2 h figures is not stated in the abstract; the protocol's benchmark section was not accessible, so those numbers are treated as paper-reported only.

## Equivalence and caveats

- **changed-output**: GRNBoost2 is not a re-implementation of GENIE3's estimator. Random forests average 1,000 bootstrapped trees with `max_features='sqrt'`; SGBM fits a boosted sequence with learning rate 0.01, 10% feature sampling and early stopping. Edge importances, the selected top-N regulons, and therefore the pruned regulons and AUCell scores all differ. Agreement should be judged by network-recovery metrics and known biology, not by matrix equality.
- The paper's speed claim is for the *pipeline*; it bundles the estimator change, parallelism and the C++/feather cisTarget rewrite, so the 10x cannot be attributed to any single step.
- The two parallel paths differ in data movement: AUCell explicitly shares one read-only ranking block (`RawArray`), whereas the GRN pool passes `ex_matrix` and `tf_matrix` inside the per-task argument tuple (`partial(...)` above). With the fork start method workers inherit those pages, but under `spawn` (an option added in 0.12.1) they are serialized per task; on very large matrices the GRN step is therefore the place where data movement must be checked, not assumed.
- Early stopping makes the tree count data-dependent; results are seed-sensitive unless `seed` is fixed, and the importance normalization divides by the realized number of trees (`to_feature_importances`).
- `arboreto`'s GENIE3 profile is "an improved implementation", not the original R GENIE3; comparing pySCENIC-genie3 to R SCENIC still changes the implementation.

## Reuse plan for Bio-Opt

- Freeze the scientific criterion **before** swapping the estimator: e.g. DREAM5 AUPRC or recovery of a known regulon set at fixed top-N. A 10x pipeline gain that loses regulons is not an optimization.
- Then test in this order:  hoist the shared predictor matrix and confirm it is read-only across tasks;  parallelize per target and measure speedup versus worker count, watching for per-task serialization of the shared matrix;  only then tune the early-stopping window and learning rate, reporting the accuracy delta alongside the time.
- Report time-to-quality, not time-to-completion: the boosted model reaches its plateau earlier, and the early-stopping window controls that trade-off directly.
