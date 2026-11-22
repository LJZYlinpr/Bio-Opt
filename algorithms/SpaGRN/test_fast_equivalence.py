"""Equivalence check: upstream bivariate Geary/Moran vs transport-only fast path."""
import json

import anndata as ad
import numpy as np
import pandas as pd


from spagrn import corexp
from spagrn.autocor import flat_weights, get_w, neighbors_and_weights

import spagrn_fast_bivariate as fast

rng = np.random.default_rng(0)
n_cells, n_genes = 300, 8
X = rng.poisson(3.0, size=(n_cells, n_genes)).astype(np.float32)
adata = ad.AnnData(X)
adata.var_names = [f"g{i}" for i in range(n_genes)]
adata.obs_names = [f"c{i}" for i in range(n_cells)]
adata.obsm["spatial"] = rng.normal(size=(n_cells, 2))
adata.layers["counts"] = X.copy()

ind, neighbors, weights_n = neighbors_and_weights(adata, latent_obsm_key="spatial", n_neighbors=5)
fw = flat_weights(adata.obs_names, ind, weights_n, n_neighbors=5)
tfs = ["g0", "g1", "g2"]
targets = ["g3", "g4", "g5", "g6", "g7"]

official_c = corexp.global_bivariate_gearys_C(
    adata.copy(), fw.copy(), tfs, targets, layer_key="counts", num_workers=2
)
fast_c = fast.geary_fast(adata.copy(), fw.copy(), tfs, targets, layer_key="counts", num_workers=2)

official_m = corexp.global_bivariate_moran_R(
    adata.copy(), fw.copy(), tfs, targets, layer_key="counts", num_workers=2
)
fast_m = fast.moran_fast(adata.copy(), fw.copy(), tfs, targets, layer_key="counts", num_workers=2)

def cmp(a, b):
    same_index = a[["TF", "target"]].equals(b[["TF", "target"]])
    max_abs = float(np.max(np.abs(a["importance"].to_numpy() - b["importance"].to_numpy())))
    exact = bool(np.array_equal(a["importance"].to_numpy(), b["importance"].to_numpy()))
    return {
        "same_tf_target_order": bool(same_index),
        "max_abs_diff": max_abs,
        "bitwise_equal": exact,
        "n_rows": [int(len(a)), int(len(b))],
    }

res = {"geary": cmp(official_c, fast_c), "moran": cmp(official_m, fast_m)}
print(json.dumps(res, indent=2))
assert all(result["bitwise_equal"] and result["same_tf_target_order"] for result in res.values())
print("EQUIVALENCE_OK")
