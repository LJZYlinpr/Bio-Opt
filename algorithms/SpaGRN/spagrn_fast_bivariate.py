"""Transport-only acceleration for SpaGRN's bivariate spatial statistics.

The upstream `global_bivariate_gearys_C` / `global_bivariate_moran_R` build a
`pool_args` list in which every single (TF, target) task carries the full dense
expression matrix (~0.5 GiB here). multiprocessing pickles that matrix once per
task, which caps throughput at a few tasks per second (measured ~4.2 tasks/s,
2,088,160 tasks -> ~5.7 days).

This module keeps the *exact* upstream worker function (`compute_pairs_c`,
which calls `bv_geary_c` / `bv_moran_r`) and the exact upstream id/weight
construction. Only the transport changes: the matrix and the flat weight arrays
are passed once per worker via the pool initializer, and tasks carry just the
two gene indices. Results are returned in the same (TF outer, target inner)
order, so the adjacency DataFrame is identical.
"""
import multiprocessing as mp

import numpy as np
import pandas as pd
import scipy

from spagrn import corexp

_G = {}
N_WORKERS = 12
CHUNK = 512


def _init(mtx, cell_x, cell_y, wij, gene_names, kind):
    _G["kind"] = kind
    _G["mtx"] = mtx
    _G["cx"] = cell_x
    _G["cy"] = cell_y
    _G["wij"] = wij
    _G["names"] = gene_names


def _work(chunk):
    fn = corexp.compute_pairs_c if _G["kind"] == "geary" else corexp.compute_pairs_m
    args = (_G["names"], _G["cx"], _G["cy"], _G["wij"], _G["mtx"])
    out = []
    for ix, iy in chunk:
        out.append(fn((args[0], ix, iy, args[1], args[2], args[3], args[4])))
    return out


def _prepare(adata, weights, tfs_in_data, select_genes, layer_key):
    gene_names = adata.var.index
    tf_ids = adata.var.index.get_indexer(tfs_in_data)
    target_ids = adata.var.index.get_indexer(select_genes)
    tmp_obs = adata.obs.copy()
    tmp_obs["id"] = np.arange(len(adata.obs))
    cell_x_id = tmp_obs.loc[weights["Cell_x"].to_list()]["id"].to_numpy()
    cell_y_id = tmp_obs.loc[weights["Cell_y"].to_list()]["id"].to_numpy()
    if layer_key:
        mtx = (
            adata.layers[layer_key].toarray()
            if scipy.sparse.issparse(adata.layers[layer_key])
            else adata.layers[layer_key]
        )
    else:
        mtx = adata.X.toarray() if scipy.sparse.issparse(adata.X) else adata.X
    wij = weights["Weight"].to_numpy()
    return gene_names, tf_ids, target_ids, cell_x_id, cell_y_id, wij, mtx


def _run(adata, weights, tfs_in_data, select_genes, layer_key, num_workers, kind="geary"):
    gene_names, tf_ids, target_ids, cx, cy, wij, mtx = _prepare(
        adata, weights, tfs_in_data, select_genes, layer_key
    )
    pairs = [(int(a), int(b)) for a in tf_ids for b in target_ids]
    chunks = [pairs[i : i + CHUNK] for i in range(0, len(pairs), CHUNK)]
    print(
        f"[fast-bivariate] kind={kind} pairs={len(pairs)} chunks={len(chunks)} workers={num_workers}",
        flush=True,
    )
    with mp.Pool(
        num_workers,
        initializer=_init,
        initargs=(mtx, cx, cy, wij, gene_names, kind),
    ) as pool:
        nested = pool.map(_work, chunks)
    rows = [r for sub in nested for r in sub]
    return pd.DataFrame(rows, columns=["TF", "target", "importance"])


def geary_fast(adata, weights, tfs_in_data, select_genes, layer_key="raw_counts", num_workers=4):
    return _run(adata, weights, tfs_in_data, select_genes, layer_key, num_workers or N_WORKERS, "geary")


def moran_fast(adata, weights, tfs_in_data, select_genes, layer_key="raw_counts", num_workers=4):
    return _run(adata, weights, tfs_in_data, select_genes, layer_key, num_workers or N_WORKERS, "moran")


def install(workers=N_WORKERS, prune_workers=6):
    global N_WORKERS
    N_WORKERS = workers
    import spagrn.regulatory_network as rn

    rn.global_bivariate_gearys_C = geary_fast
    rn.global_bivariate_moran_R = moran_fast

    orig_receptor_auc = rn.InferNetwork.receptor_auc

    def receptor_auc(self, *args, **kwargs):
        kwargs.setdefault("num_workers", workers)
        return orig_receptor_auc(self, *args, **kwargs)

    rn.InferNetwork.receptor_auc = receptor_auc

    orig_prune = rn.InferNetwork.prune_modules

    def prune_modules(self, modules, dbs, motif_anno_fn, num_workers, **kwargs):
        n = min(int(num_workers or prune_workers), prune_workers)
        print(f"[resource] prune_modules workers={n}", flush=True)
        return orig_prune(self, modules, dbs, motif_anno_fn, n, **kwargs)

    rn.InferNetwork.prune_modules = prune_modules
    print(
        f"[fast-bivariate] installed (bivariate_workers={workers}, prune_workers={prune_workers})",
        flush=True,
    )
