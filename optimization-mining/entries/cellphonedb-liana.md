---
slug: cellphonedb-liana
method: CellPhoneDB
omicverse_binding: "omicverse.single.run_liana (imports the liana package, default method='rank_aggregate'), omicverse.single.cellphonedb_v5 / run_cellphonedb_v5, and ov.single.run_ccc(method='liana'|'cellphonedb'); exported in omicverse v2.4.3"
original:
  paper: "Efremova et al., Nat Protoc 15, 1484-1506 (2020), doi:10.1038/s41596-020-0292-x"
  repo: "https://github.com/ventolab/CellphoneDB"
  commit: "dc8abd15b24e1d48c8862d54e8122603630a68ed"
optimized:
  paper: "Dimitrov et al., Nat Cell Biol 26, 1613-1622 (2024), doi:10.1038/s41556-024-01469-w"
  repo: "https://github.com/saezlab/liana-py"
  commit: "6bc8da9f7f0207d809b0702bb68789efb6419a73"
optimization_type: system
mechanism: "Reimplement each LR method, including CellPhoneDB, as a shared set of sparse statistics scored by compiled parallel permutation kernels, then aggregate methods by rank instead of running each tool separately."
reported_speedup: "Not reported as an explicit CellPhoneDB-to-LIANA+ ratio. LIANA+ Extended Data Fig. 1 benchmarks 'scalable implementations' of the dissociated methods on simulated data from 1,000 to 100,000 cells, reporting time in log2 seconds and peak RAM; the text claims scalable implementations but gives no numeric fold change."
equivalence: changed-output
optimization_families: [O4, O6, O9, O10]
sources:
  - "https://doi.org/10.1038/s41596-020-0292-x"
  - "https://github.com/ventolab/CellphoneDB/blob/master/README.md"
  - "https://github.com/ventolab/CellphoneDB/blob/master/cellphonedb/src/core/methods/cpdb_statistical_analysis_helper.py"
  - "https://github.com/ventolab/CellphoneDB/blob/master/cellphonedb/src/core/methods/cpdb_statistical_analysis_method.py"
  - "https://doi.org/10.1038/s41556-024-01469-w"
  - "https://doi.org/10.1038/s41467-022-30755-0"
  - "https://github.com/saezlab/liana-py/blob/main/src/liana/_core/_pipe_utils/_get_mean_perms.py"
  - "https://github.com/saezlab/liana-py/blob/main/src/liana/method/sc/_cellphonedb.py"
  - "https://github.com/saezlab/liana-py/blob/main/src/liana/method/sc/_liana_pipe.py"
verified: partial
---

## Baseline bottleneck

CellPhoneDB's statistical analysis shuffles cell-type labels and recomputes, for every permutation, the per-cluster mean of every gene (detection percentages are skipped in this mode), then scores every interaction and cluster pair against the observed values. Defaults are `iterations=1000`, `threads=4`, `threshold=0.1` (`cpdb_statistical_analysis_method.py` L15–17; README L148–149, L163). The p-value floor is 1/1000. The pinned commit is the current CellPhoneDB release line; the 2020 paper documents the same shuffle-and-recompute test, while the current README carries the defaults used here (https://github.com/ventolab/CellphoneDB/blob/master/README.md).

The hot path rebuilds a dense cluster-mean table each iteration inside Python/pandas:

```python
def shuffled_analysis(iterations, meta, counts, ...):
    ...
    with Pool(processes=threads) as pool:
        for result in pool.imap(statistical_analysis_thread, range(iterations)):
            results.append(result)

def _statistical_analysis(..., iteration_number):
    shuffled_meta = shuffle_meta(meta)
    shuffled_clusters = build_clusters(shuffled_meta, counts,
                                       complex_to_protein_ids, skip_percent=True)
    shuffled_mean_analysis = mean_analysis(interactions, shuffled_clusters,
                                           cluster_combinations, separator)
    return shuffled_greater_than_real(shuffled_mean_analysis, real_mean_analysis)
```

`build_clusters` calls `npg.aggregate(codes, counts.values, func='mean', axis=1)` over the full genes x cells matrix and wraps the result in a `pandas.DataFrame` (`cpdb_statistical_analysis_helper.py` L101–169). So the cost is `iterations * (full-matrix aggregation + DataFrame construction + all cluster-pair scoring)`, with only process-level parallelism and a default of 4 workers.

## Optimization trajectory

1. **Assemble the ligand/receptor statistics once .** LIANA builds the shared per-entity table (ligand/receptor means, detection proportions, complex assignments) once per run, and every method re-scores that table rather than re-reading the matrix (`_liana_pipe.liana_pipe_consensus`, L340–450). A nine-method consensus therefore pays data preparation once. The permutation null cube is still regenerated per permutation-based method; only the per-entity statistics are shared.
2. **Move the permutation aggregation into a compiled sparse kernel .** `_perm_group_sums` walks the CSR buffers directly, adds each cell's non-zeros to the label it carries in the permuted assignment, and returns the `(n_perms, n_labels, n_genes)` null cube (`_get_mean_perms.py` L177–200, https://github.com/saezlab/liana-py/blob/main/src/liana/_core/_pipe_utils/_get_mean_perms.py). No permuted copy of the matrix is materialized and no pandas object is built per permutation; the same file has a trimean variant that buckets non-zeros per gene and sorts each gene's column (`_perm_group_trimeans`).
3. **Parallelize over permutations .** The kernels are `@nb.njit(parallel=True)` and loop with `nb.prange`; `n_jobs` maps to the numba thread count, and permutation order is independent of `n_jobs`, so p-values depend only on the seed (`_numba_threads`, `_generate_perms_cube`).
4. **Chunk to bound memory .** Permutations are drawn and aggregated a block at a time (`_chunk_permutations`, `_MAX_PERM_INDEX_ELEMENTS`, `_MAX_PERM_STAT_ELEMENTS`), and p-values are accumulated blockwise in `_perm_pvals`, so peak memory does not scale with `n_perms` (L364–410).
5. **Vectorize scoring and p-value counting .** Only the ligand/receptor columns at the interaction's source/target labels are gathered from the null cube, the method's score is applied elementwise to the block, and exceedances are counted with a tie tolerance (`_perm_pvals`). This avoids materializing the full `(n_perms, n_interactions)` matrix at once (only bounded blocks) and avoids Python loops over interactions.
6. **Aggregate methods instead of running tools separately (workflow).** The default `rank_aggregate` runs nine methods on the same statistics and combines them by robust rank aggregation, with CellPhoneDB as one reimplemented method (`_cellphonedb.py`: mean of ligand and receptor means, zero when either is absent; `permute=True`).
7. **Extend to spatial and multi-condition data.** LIANA+ adds spatially weighted local metrics, multi-view factor analysis and multi-sample designs in the same framework (paper, Figs. 1–2).

```python
@nb.njit(parallel=True, cache=True)
def _perm_group_sums(data, indices, indptr, perm_indices,
                     label_of_row, n_labels, n_genes):
    n_perms, n_obs = perm_indices.shape
    sums = np.zeros((n_perms, n_labels, n_genes), dtype=np.float64)
    for p in prange(n_perms):
        for j in range(n_obs):
            label = label_of_row[j]
            row = perm_indices[p, j]
            for k in range(indptr[row], indptr[row + 1]):
                sums[p, label, indices[k]] += data[k]
    return sums
```

Why this reduces work: the permutation count is unchanged (1,000 by default in both), but each permutation becomes one sparse pass in compiled code instead of a full dense `npg.aggregate` plus DataFrame construction, and scoring touches only the O(n_interactions) entries. The asymptotic `O(n_perms * nnz)` remains; the gain is the constant factor, the removal of pandas overhead, and bounded peak memory.

## Evidence

- LIANA+ Extended Data Fig. 1: "LIANA+ provides scalable implementations of cell-cell communication methods. Efficiency benchmark of A) dissociated methods and B) spatially weighted local metrics implemented in LIANA+. Time in log2 of seconds and peak total RAM (memory) usage are shown over a number of simulated datasets with observations (cells) ranging from 1000 to 100,000" (paper, Extended Data Fig. 1 caption). No numeric fold change is stated in the text.
- LIANA (2022) Methods: "CellPhoneDB's algorithm was re-implemented in LIANA and used throughout this manuscript with 1000 permutations. Identical to the original implementation, cluster labels were reshuffled and an one-sided empirical p-value was calculated..." (doi:10.1038/s41467-022-30755-0). This establishes the reimplementation lineage and the same permutation count.
- Code: compiled CSR permutation kernels, blockwise p-values and shared statistics as cited above (`_get_mean_perms.py`, `_liana_pipe.py`, `_cellphonedb.py`).
- Defaults: `V.n_perms = 1000`, `V.expr_prop = 0.05` (`liana/_core/_constants.py`), while CellPhoneDB uses `threshold=0.1` and `iterations=1000`.

## Equivalence and caveats

`changed-output`. The default LIANA+ entry point (`rank_aggregate`) returns a consensus rank across nine methods, not CellPhoneDB's per-interaction means and p-values; its columns and significance calls are different by construction. The per-method path (`method='cellphonedb'`) is the closest analogue: it uses the same score (mean of ligand and receptor means, zero if either is absent) and a one-sided label-permutation p-value, but it is a reimplementation, not a call into CellPhoneDB, and defaults differ (`expr_prop=0.05` vs `threshold=0.1`), so the tested interaction set and p-values can differ. Complex handling and resource versions also differ (LIANA uses its packaged consensus/OmniPath resources by default, CellPhoneDB ships its own database). Multi-testing correction and rank aggregation further change which interactions are reported.

Do not treat a LIANA+ consensus score as a CellPhoneDB p-value. For a like-for-like comparison, set `expr_prop=0.1`, the same resource/complex table, the same `n_perms` and seed, and use `method='cellphonedb'`.

## omicverse binding

`omicverse.single.run_liana` imports `liana` and calls `li.mt.<method>` with `method='rank_aggregate'` by default, storing `adata.uns['liana_res']`; it requires the optional `liana` package (`omicverse/single/_liana.py` L456–525). `omicverse.single.cellphonedb_v5` / `run_cellphonedb_v5` call `cpdb_statistical_analysis_method.call` and are the original-tool path, not LIANA (omicverse defaults: `iterations=1000`, `threshold=0.1`, `threads=10`; CellPhoneDB's own default is `threads=4`). `ov.single.run_ccc(method='liana'|'cellphonedb')` dispatches to these two backends (`omicverse/single/_ccc.py` L37–59, L293). All are exported in v2.4.3 (`omicverse/single/__init__.py` L65–84, L370–391).

## Reuse plan for Bio-Opt

Patterns: ** reuse shared state**, ** vectorize/compile**, ** chunk**, ** parallelize**.

Test first: on one dataset with the same resource and `expr_prop=0.1`, compare LIANA's `cellphonedb` output with CellPhoneDB v5 (means, p-values, and the significant set at p<0.05) and measure runtime versus `n_jobs`. Then profile the permutation kernel to confirm that the residual time is the sparse pass rather than scoring or DataFrame assembly. Keep the consensus output out of any equivalence claim until the per-method path is validated.
