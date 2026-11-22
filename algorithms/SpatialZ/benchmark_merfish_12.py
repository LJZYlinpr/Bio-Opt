import argparse
import importlib.util, inspect, json, sys, time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
import torch

parser = argparse.ArgumentParser(description="SpatialZ exact-output 12-layer MERFISH comparison")
parser.add_argument("--input-dir", type=Path, required=True)
parser.add_argument("--out", type=Path, required=True)
args = parser.parse_args()
args.out.mkdir(parents=True, exist_ok=False)
root = Path(__file__).resolve().parent
sys.path.insert(0, str(root / "optimized"))
from cpu import measure, synthesize_cpp
from full_pair import assign_labels

NATIVE = root / "original/SpatialZ.py"

def dense(x):
    return x.toarray() if sp.issparse(x) else np.asarray(x)

def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m

native = load_module("spatialz_native", NATIVE)
fast = load_module("geometry_only_spatialz", NATIVE)
source = NATIVE.read_text()
marker = "    Coor_final = Coor_torch.detach().cpu().numpy()"
exec(compile(source.replace(marker, marker + "\n    return Coor_final", 1), str(NATIVE), "exec"), fast.__dict__)
src = inspect.getsource(synthesize_cpp)
scope = dict(synthesize_cpp.__globals__); scope["demo_beta"] = 1
exec(src.replace("100 * similarities", "demo_beta * similarities", 1), scope)
synthesize_cpp = scope["synthesize_cpp"]

directory = args.input_dir
paths = sorted(directory.glob("merfish_*_paste.h5ad"))
if len(paths) != 5:
    raise ValueError("Expected the five official MERFISH H5AD files")
data = ad.concat([ad.read_h5ad(f) for f in paths], join="outer", merge="same")
keys = ["0.04", "0.09", "0.14", "0.19", "0.24"]
frames = [data[data.obs.slice_id.astype(str) == k].copy() for k in keys]
for f in frames:
    f.obsm["spatial"] = np.asarray(f.obsm["spatial"])[:, :2]
    f.uns.clear()
    for key in list(f.obsp):
        del f.obsp[key]
config = {"lr": 1e-5, "num_projections": 80, "Beta": 1}

def prep(f):
    g = f.copy()
    if sp.issparse(g.X):
        g.X = g.X.toarray()
    return g

pairs = [(0, 1), (1, 2), (2, 3), (3, 4)]
res = {"contract": "np.random.seed(42) before each stock call; synth_cpp(seed=42)", "pairs": []}
t_stock_all = t_ours_all = 0.0
for (l, r) in pairs:
    a, b = frames[l], frames[r]
    id1, id2 = keys[l], keys[r]
    t0 = time.time()
    stock = {}
    for j in (1, 2, 3):
        alpha = 1 - j / 4
        np.random.seed(42)
        sim = native.Generate_spatialz(
            a.copy(), b.copy(), adata1_id=id1, adata2_id=id2, alpha=alpha, device="cuda:0",
            n_cell=None, n_mag=1.0, lr=1e-5, nb_iter_max=1000, seed=42, num_projections=80,
            cell_type_key="cell_class", syn_mode="default", k_sam=50, micro_env_key="mender",
            Beta=1, add_obs_list=None, verbose=False,
        )
        stock[j] = (dense(sim.X), np.asarray(sim.obsm["spatial"]), sim.obs["cell_class"].astype(str).to_numpy())
    t_stock = time.time() - t0

    t1 = time.time()
    qs = []
    for j in (1, 2, 3):
        alpha = 1 - j / 4
        kw = dict(alpha=alpha, n_cell=None, nb_iter_max=1000, seed=42, n_mag=1.0,
                  syn_mode="default", k_sam=50, device="cuda:0", cell_type_key="cell_class", **config)
        xy = fast.Generate_spatialz(a.copy(), b.copy(), **kw)
        labels = assign_labels(a, b, xy)
        q = ad.AnnData(np.zeros((len(xy), a.n_vars), dtype=np.float32),
                       obs=pd.DataFrame({"cell_class": pd.Categorical(labels)}), var=a.var.copy())
        q.obsm["spatial"] = xy
        qs.append((j, q, labels))
    allf = [prep(a), prep(b)] + [prep(q) for _, q, _ in qs]
    d = measure(allf, 32)
    for f, m in zip(allf, d):
        f.obsm["mender"] = m
    ours = {}
    for i, (j, q, labels) in enumerate(qs):
        arr, _ = synthesize_cpp([allf[0], allf[1], allf[2 + i]], seed=42, threads=4)
        ours[j] = (np.asarray(arr), np.asarray(q.obsm["spatial"]), np.asarray(labels))
    t_ours = time.time() - t1

    layers = []
    for j in (1, 2, 3):
        se, sxy, slb = stock[j]
        oe, oxy, olb = ours[j]
        oe = oe.astype(se.dtype) if oe.dtype != se.dtype else oe
        layers.append({"j": j,
                       "expr_exact": bool(np.array_equal(oe, se)),
                       "expr_diff": int(np.count_nonzero(oe != se)),
                       "xy_exact": bool(np.array_equal(oxy, sxy)),
                       "labels_exact": bool(np.array_equal(olb, slb))})
    res["pairs"].append({"pair": [id1, id2], "stock_s": round(t_stock, 3), "ours_s": round(t_ours, 3),
                         "layers": layers})
    t_stock_all += t_stock
    t_ours_all += t_ours
    print(json.dumps(res["pairs"][-1]), flush=True)
    p = args.out
    (p / "seeded_full_RESULT.json").write_text(json.dumps(res, indent=2) + "\n")

res["stock_total_s"] = round(t_stock_all, 3)
res["ours_total_s"] = round(t_ours_all, 3)
res["speedup"] = round(t_stock_all / t_ours_all, 3)
res["all_layers_exact"] = all(x["expr_exact"] and x["xy_exact"] and x["labels_exact"]
                              for p in res["pairs"] for x in p["layers"])
p = args.out
(p / "seeded_full_RESULT.json").write_text(json.dumps(res, indent=2) + "\n")
print(json.dumps({k: res[k] for k in ("stock_total_s", "ours_total_s", "speedup", "all_layers_exact")}), flush=True)

if not res["all_layers_exact"]:
    raise SystemExit("SpatialZ output equivalence failed")
