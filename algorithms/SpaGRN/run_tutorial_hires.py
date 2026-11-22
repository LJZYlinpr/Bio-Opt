"""Official SpaGRN high-resolution tutorial (stereo_seq_mouse_brain_hi-res.ipynb).

Reproduces the notebook cells in order. Only these host/input adaptations:
  * layers='counts' (the public h5ad ships raw counts in layers['counts'],
    the notebook text says 'raw_counts')
  * num_workers=12 (12-CPU quota; notebook uses 20)
  * explicit output_dir, and matplotlib is headless (figures saved, not shown)
"""
import json
import os
import sys
import resource
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from spagrn.regulatory_network import InferNetwork as irn

R = Path(os.environ["R_ROOT"])
DATA = R / "data/SpaGRN"
WORK = R / "work/SpaGRN"
OUT = WORK / "exp_output"
OUT.mkdir(parents=True, exist_ok=True)
os.chdir(WORK)

N_WORKERS = 12
FAST_BIVARIATE = os.environ.get("SPAGRN_FAST_BIVARIATE") == "1"
if FAST_BIVARIATE:
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import spagrn_fast_bivariate as _fast

    _fast.install(N_WORKERS, int(os.environ.get("SPAGRN_PRUNE_WORKERS", "6")))
marks = {}


def mark(name, t0):
    marks[name] = round(time.time() - t0, 2)
    print(f"[stage] {name} {marks[name]:.1f}s", flush=True)


def main():
    t_all = time.time()

    # --- cell 4: load data -------------------------------------------------------
    t0 = time.time()
    fn = str(DATA / "Mouse_brain_cell_bin.h5ad")
    data = irn.read_file(fn)
    print(data, flush=True)
    shape_loaded = list(data.shape)
    mark("load", t0)

    # --- cell 6: preprocess ------------------------------------------------------
    t0 = time.time()
    data = irn.preprocess(data, min_genes=1000, min_cells=300, min_counts=10, max_gene_num=4000)
    print(data, flush=True)
    shape_pre = list(data.shape)
    mark("preprocess", t0)

    # --- cell 8: PCA -------------------------------------------------------------
    t0 = time.time()
    sc.tl.pca(data)
    sc.pl.pca(data, color="annotation", show=False, save="_tutorial_hires.png")
    mark("pca", t0)

    # --- cell 10: TF list and databases -----------------------------------------
    tfs_fn = str(DATA / "mouse_tfs.txt")
    database_fn = str(DATA / "mm10_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather")
    motif_anno_fn = str(DATA / "motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl")

    # --- cell 12-14: object and parameters --------------------------------------
    grn = irn(data, project_name="mouse_brain")
    grn.add_params({"prune_auc_threshold": 0.05, "rank_threshold": 9000, "auc_threshold": 0.05})

    # --- cell 16: ligand-receptor ------------------------------------------------
    niche_human = pd.read_csv(DATA / "resource/lr_network_human.csv")
    niche_mouse = pd.read_csv(DATA / "resource/lr_network_mouse.csv")
    niches = pd.concat([niche_mouse, niche_human])
    print(niches, flush=True)

    # --- cell 18: construct the network -----------------------------------------
    t0 = time.time()
    adata = grn.infer(
        database_fn,
        motif_anno_fn,
        tfs_fn,
        niche_df=niches,
        gene_list=None,
        num_workers=N_WORKERS,
        cache=True,
        output_dir=str(OUT),
        save_tmp=True,
        layers="counts",
        latent_obsm_key="spatial",
        model="bernoulli",
        n_neighbors=10,
        methods=["FDR_I", "FDR_C", "FDR_G"],
        operation="intersection",
        mode="geary",
        cluster_label="annotation",
    )
    mark("infer", t0)

    # --- cells 20-26: reload written result and validate ------------------------
    t0 = time.time()
    res_fn = OUT / "mouse_brain_spagrn.h5ad"
    adata = sc.read_h5ad(res_fn)
    print(adata, flush=True)

    regulons = adata.uns.get("regulon_dict")
    receptors = adata.uns.get("receptor_dict")
    isr = adata.obsm.get("isr")
    auc = adata.obsm.get("auc_mtx")

    summary = {
        "output": str(res_fn),
        "shape_loaded": shape_loaded,
        "shape_after_preprocess": list(adata.shape),
        "n_obs": int(adata.n_obs),
        "n_vars": int(adata.n_vars),
        "uns_keys": list(adata.uns.keys()),
        "obsm_keys": list(adata.obsm_keys()),
        "n_regulons": int(len(regulons)) if regulons is not None else 0,
        "regulons": sorted(map(str, regulons)) if regulons is not None else [],
        "regulon_sizes": ({str(k): int(len(v)) for k, v in regulons.items()} if regulons is not None else {}),
        "n_receptor_regulons": int(len(receptors)) if receptors is not None else 0,
        "receptor_dict_keys": sorted(map(str, receptors)) if receptors is not None else [],
        "isr_shape": list(isr.shape) if isr is not None else None,
        "auc_mtx_shape": list(auc.shape) if auc is not None else None,
        "timings_s": marks,
        "total_s": round(time.time() - t_all, 2),
        "peak_rss_kb": int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        "child_peak_rss_kb": int(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss),
        "n_workers": N_WORKERS,
        "layer_key_used": "counts",
        "bivariate_transport": "fast_shared_matrix" if FAST_BIVARIATE else "upstream_pool_args",
    }
    (WORK / "run_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2), flush=True)
    mark("validate", t0)
    print("SPAGRN_TUTORIAL_HIRES_DONE", flush=True)


if __name__ == "__main__":
    main()
