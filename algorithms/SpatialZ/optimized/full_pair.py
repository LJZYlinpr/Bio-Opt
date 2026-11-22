"""Run one stage of an ABCA-1 adjacent-depth pair, with per-plane checkpoints."""
import argparse
import ctypes
import gc
import json
from pathlib import Path
import time

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from sklearn.neighbors import NearestNeighbors

from cpu import atomic_json, atomic_npz
from layers import target_layers

ROOT = Path(__file__).resolve().parent / "full"


def process_rss_mib():
    for line in Path("/proc/self/status").read_text().splitlines():
        if line.startswith("VmRSS:"):
            return round(int(line.split()[1]) / 1024, 2)
    raise RuntimeError("Unable to read current RSS")


def geometry(pair):
    import h5py
    import torch
    from geometry import CONFIG, NATIVE, load_module

    torch.set_num_threads(1)
    if torch.ones(2, device="cuda:0").sum().item() != 2:
        raise RuntimeError("CUDA validation failed")
    native = load_module("spatialz_full_native", NATIVE)
    source = NATIVE.read_text()
    line = "    Coor_final = Coor_torch.detach().cpu().numpy()"
    assert source.count(line) == 1
    source = source.replace(line, line + "\n    return Coor_final", 1)
    exec(compile(source, str(NATIVE), "exec"), native.__dict__)
    a, b = pair["left"], pair["right"]
    frames = []
    for endpoint in (a, b):
        with h5py.File(ROOT / "sections" / f'{endpoint["id"]}.h5ad') as f:
            xy = f["obsm/spatial"][:]
        frame = ad.AnnData(csr_matrix((len(xy), 0), dtype=np.float32),
                           obs=pd.DataFrame({"cell_class": np.zeros(len(xy), dtype=np.int8)}))
        frame.obsm["spatial"] = xy
        frames.append(frame)
    for key, z, _ in target_layers(pair):
        path = ROOT / "geometry" / pair["id"] / f"{key}.npz"
        if path.exists():
            continue
        torch.cuda.reset_peak_memory_stats()
        started = time.monotonic()
        alpha = 1 - (z - a["z_um"]) / (b["z_um"] - a["z_um"])
        xy = native.Generate_spatialz(
            frames[0].copy(), frames[1].copy(), adata1_id=a["id"], adata2_id=b["id"],
            alpha=alpha, cell_type_key="cell_class", device="cuda:0", verbose=True, **CONFIG,
        )
        if xy.ndim != 2 or xy.shape[1] != 2 or not np.isfinite(xy).all():
            raise ValueError(f"Invalid geometry for z={z}")
        path.parent.mkdir(parents=True, exist_ok=True)
        atomic_npz(path, xy=xy, z_um=z, seconds=time.monotonic() - started)
        n_cells = len(xy)
        del xy
        torch.cuda.empty_cache()
        print("GEOMETRY", pair["id"], z, n_cells, time.monotonic() - started,
              "GPU_MEMORY_MIB", {"peak_allocated": torch.cuda.max_memory_allocated() / 2**20,
                                 "after_allocated": torch.cuda.memory_allocated() / 2**20,
                                 "after_reserved": torch.cuda.memory_reserved() / 2**20}, flush=True)


def assign_labels(a, b, xy):
    labels = []
    distance = []
    for frame in (a, b):
        model = NearestNeighbors(n_neighbors=1).fit(frame.obsm["spatial"])
        ds, index = model.kneighbors(xy)
        labels.append(frame.obs["cell_class"].to_numpy(dtype=str)[index[:, 0]])
        distance.append(ds[:, 0])
    # The author uses the inverse of distance plus 0.1 and sorts group names for an exact tie.
    w_left, w_right = 1 / (distance[0] + 0.1), 1 / (distance[1] + 0.1)
    left = (w_left > w_right) | ((w_left == w_right) & (labels[0] <= labels[1]))
    return np.where(left, labels[0], labels[1])


def wait_for_geometry(pair_id, key):
    path = ROOT / "geometry" / pair_id / f"{key}.npz"
    state_path = ROOT / "queue_state" / f"gpu-{pair_id}.json"
    while not path.is_file():
        if state_path.is_file() and json.loads(state_path.read_text()).get("state") == "FAILED":
            raise RuntimeError(f"GPU pair {pair_id} failed before geometry for layer {key}")
        time.sleep(1)
    return path


def expression_plane(a, b, genes, pair_id, z, seed, path, out):
    from cpu import measure, synthesize_cpp

    with np.load(path) as data:
        xy = data["xy"]
        if "z_um" in data and data["z_um"].item() != z:
            raise ValueError(f"Geometry layer depth in {path} does not match {z}")
    q = ad.AnnData(csr_matrix((len(xy), len(genes)), dtype=np.float32),
                   obs=pd.DataFrame({"cell_class": pd.Categorical(assign_labels(a, b, xy))}),
                   var=pd.DataFrame(index=genes))
    q.obsm["spatial"] = xy
    frames = [a, b, q]
    features = measure(frames, 32)
    for frame, values in zip(frames, features):
        frame.obsm["mender"] = values
    expr, stats = synthesize_cpp(frames, seed=seed, threads=4)
    assert expr.shape == q.shape and np.isfinite(expr).all() and (expr >= 0).all()
    q.X = csr_matrix(expr)
    q.uns["z_um"] = z
    q.uns["expression_scale"] = "raw_counts"
    q.uns["source_pair"] = pair_id
    q.uns["rng_seed"] = seed
    for frame in frames:
        del frame.obsm["mender"]
    out.parent.mkdir(parents=True, exist_ok=True)
    temp = out.with_suffix(".tmp.h5ad")
    q.write_h5ad(temp, compression="lzf")
    temp.replace(out)
    stats["rss_before_release_mib"] = process_rss_mib()
    return q.shape, stats


def expression(pair, genes):
    a, b = [ad.read_h5ad(ROOT / "sections" / f'{p["id"]}.h5ad') for p in (pair["left"], pair["right"])]
    assert a.var_names.tolist() == b.var_names.tolist() == genes
    # The sample kernel expects float32 dense counts. MENDER receives the native
    # AnnData with the same count matrix, labels and spatial coordinates.
    for frame in (a, b):
        frame.X = frame.X.toarray().astype(np.float32, copy=False)
    trim = ctypes.CDLL(None).malloc_trim
    trim.argtypes = [ctypes.c_size_t]
    trim.restype = ctypes.c_int
    for key, z, seed in target_layers(pair):
        out = ROOT / "predictions" / pair["id"] / f"{key}.h5ad"
        if out.exists():
            continue
        path = wait_for_geometry(pair["id"], key)
        started = time.monotonic()
        shape, stats = expression_plane(a, b, genes, pair["id"], z, seed, path, out)
        gc.collect()
        trim(0)
        stats["rss_after_release_mib"] = process_rss_mib()
        print("EXPRESSION", pair["id"], z, shape, time.monotonic() - started, stats, flush=True)


def main():
    global ROOT
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["gpu", "cpu"])
    parser.add_argument("id")
    parser.add_argument("--root", type=Path, default=ROOT)
    args = parser.parse_args()
    ROOT = args.root
    data = json.loads((ROOT / "QUEUE.json").read_text())
    pair = next(p for p in data["pairs"] if p["id"] == args.id)
    started = time.monotonic()
    if args.stage == "gpu":
        geometry(pair)
    else:
        expression(pair, data["gene_intersection"])
    atomic_json(ROOT / f"{args.stage}_done" / f"{pair['id']}.json",
                {"pair": pair["id"], "targets": pair["targets"], "seconds": time.monotonic() - started})


if __name__ == "__main__":
    main()
