"""CPU continuation and controlled comparisons for the SpatialZ pilot."""
import argparse
import contextlib
import ctypes
import io
import json
import os
from pathlib import Path
import resource
import time

import anndata as ad
import MENDER
import numpy as np
from scipy.spatial import cKDTree
from sklearn.preprocessing import normalize

ROOT = Path(__file__).resolve().parent
NATIVE = ROOT.parent / "original/SpatialZ.py"


def atomic_json(path, data):
    temp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    temp.write_text(json.dumps(data, indent=2) + "\n")
    temp.replace(path)


def atomic_npz(path, **data):
    temp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with temp.open("wb") as stream:
        np.savez(stream, **data)
    temp.replace(path)


def load_frames(target):
    folder = ROOT / "geometry" / str(target)
    assert (folder / "geometry.json").is_file(), folder
    return [ad.read_h5ad(folder / f"{name}.h5ad") for name in ("left", "right", "query")]


def measure(frames, workers):
    # Preserve native category order and transductive descriptor computation.
    for index, frame in enumerate(frames):
        frame.obs["temp_slice_id"] = f"temp_adata{index + 1}"
    combined = frames[0].concatenate(frames[1:])
    for key in ("temp_slice_id", "cell_class"):
        combined.obs[key] = combined.obs[key].astype("category")
    model = MENDER.MENDER(combined, batch_obs="temp_slice_id", ct_obs="cell_class", random_seed=0)
    model.prepare()
    model.set_MENDER_para(n_scales=6, nn_mode="radius", nn_para=15)
    with contextlib.redirect_stdout(io.StringIO()):
        model.run_representation_mp(workers)
    combined.obsm["mender"] = model.adata_MENDER.X
    result = []
    for index, frame in enumerate(frames):
        data = combined[combined.obs["temp_slice_id"] == f"temp_adata{index + 1}"]
        result.append(np.array(data.obsm["mender"], copy=True))
        del frame.obs["temp_slice_id"]
    return result


def synthesize_native(frames, seed):
    import importlib.util
    spec = importlib.util.spec_from_file_location("spatialz_native", NATIVE)
    native = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(native)
    np.random.seed(seed)
    # The author prints one warning per row with no matching candidate.
    with contextlib.redirect_stdout(io.StringIO()):
        return native.synthesize_gene_expression(
            frames[:2], frames[2], "cell_class", 50, "mender", 100,
        ).X


def synthesize_cpp(frames, seed, threads):
    a, b, q = frames
    ref_types = np.concatenate([f.obs["cell_class"].to_numpy() for f in (a, b)])
    query_types = q.obs["cell_class"].to_numpy()
    ref_xy = np.vstack([f.obsm["spatial"] for f in (a, b)])
    tree = cKDTree(ref_xy)
    ref_desc = normalize(np.vstack([f.obsm["mender"] for f in (a, b)]))
    query_desc = normalize(q.obsm["mender"])
    width = 50
    candidates = np.zeros((q.n_obs, width), dtype=np.int64)
    cumulative = np.zeros((q.n_obs, width), dtype=np.float64)
    lengths = np.zeros(q.n_obs, dtype=np.int64)
    preparing = time.monotonic()

    for cell_type in np.unique(query_types):
        rows = np.flatnonzero(query_types == cell_type)
        k = min(width, int(np.count_nonzero(ref_types == cell_type)))
        if not k:
            continue
        _, indices = tree.query(q.obsm["spatial"][rows], k=k)
        indices = np.asarray(indices).reshape(len(rows), k)
        # Author queries the global tree first, then filters its candidates.
        match = ref_types[indices] == cell_type
        local_row, local_col = np.nonzero(match)
        if not len(local_row):
            continue
        donors = indices[local_row, local_col]
        owning_rows = rows[local_row]
        similarities = np.einsum("ij,ij->i", ref_desc[donors], query_desc[owning_rows])
        weights = np.exp(100 * similarities)
        offsets = np.column_stack((np.zeros(len(rows), dtype=np.int64), match.cumsum(1)))
        ranks = offsets[local_row, local_col]
        lengths[rows] = offsets[:, -1]
        candidates[owning_rows, ranks] = donors
        weight_table = np.zeros((len(rows), k), dtype=weights.dtype)
        weight_table[local_row, ranks] = weights
        totals = weight_table.sum(axis=1, keepdims=True)
        np.divide(weight_table, totals, out=weight_table, where=totals > 0)
        bounds = weight_table.cumsum(axis=1, dtype=np.float64)
        np.divide(bounds, bounds[:, -1:], out=bounds, where=bounds[:, -1:] > 0)
        cumulative[rows, :k] = bounds

    if not np.isfinite(cumulative).all():
        raise ValueError("SpatialZ candidate probabilities are nonfinite")
    probabilities_ready = time.monotonic()
    # Generate the same MT19937 sequence in native query/gene order.
    draw_rows = np.maximum(np.cumsum(lengths > 0) - 1, 0)
    random = np.random.RandomState(seed)
    uniform = random.random_sample((int(np.count_nonzero(lengths)), q.n_vars))
    reference = np.ascontiguousarray(np.vstack((a.X, b.X)), dtype=np.float32)
    output = np.zeros(q.shape, dtype=np.float32)
    integers = np.ctypeslib.ndpointer(np.int64, flags="C_CONTIGUOUS")
    doubles = np.ctypeslib.ndpointer(np.float64, flags="C_CONTIGUOUS")
    floats = np.ctypeslib.ndpointer(np.float32, flags="C_CONTIGUOUS")
    kernel = ctypes.CDLL(str(ROOT / "sample.so")).sample_expression
    kernel.argtypes = [floats, integers, doubles, integers, integers, doubles, floats,
                       ctypes.c_int64, ctypes.c_int64, ctypes.c_int64, ctypes.c_int]
    kernel.restype = None
    kernel_start = time.monotonic()
    kernel(reference, candidates, cumulative, lengths, draw_rows, uniform, output,
           q.n_obs, q.n_vars, width, threads)
    metrics = {
        "candidate_seconds": probabilities_ready - preparing,
        "cpp_kernel_seconds": time.monotonic() - kernel_start,
        "rows_without_candidates": int(np.count_nonzero(lengths == 0)),
    }
    return output, metrics


def run(args):
    started = time.monotonic()
    frames = load_frames(args.target)
    out = ROOT / "pilot" / args.tag / str(args.target)
    out.mkdir(parents=True, exist_ok=True)
    descriptor = ROOT / "pilot" / args.descriptor_tag / str(args.target) / "descriptor.npz"
    info = {"mode": args.mode, "workers": args.workers, "threads": args.threads,
            "seed": args.seed, "target": args.target, "node": os.uname().nodename,
            "job": os.environ.get("SLURM_JOB_ID"), "step": os.environ.get("SLURM_STEP_ID")}
    if args.mode in {"mender", "full"}:
        t = time.monotonic()
        arrays = measure(frames, args.workers)
        info["mender_seconds"] = time.monotonic() - t
        info["descriptor_dtype"] = [str(a.dtype) for a in arrays]
        if args.mode == "mender":
            atomic_npz(out / "descriptor.npz", **{str(i): a for i, a in enumerate(arrays)})
    else:
        with np.load(descriptor) as data:
            arrays = [data[str(i)] for i in range(3)]
    if args.mode != "mender":
        for frame, features in zip(frames, arrays):
            frame.obsm["mender"] = features
        t = time.monotonic()
        q = frames[-1]
        if args.mode == "native":
            expr = synthesize_native(frames, args.seed)
        else:
            expr, stage = synthesize_cpp(frames, args.seed, args.threads)
            info.update(stage)
        info["synthesis_seconds"] = time.monotonic() - t
        assert expr.shape == q.shape and np.isfinite(expr).all() and np.min(expr) >= 0
        atomic_npz(out / "prediction.npz", expr=expr, xy=q.obsm["spatial"],
                   labels=q.obs["cell_class"].to_numpy(dtype=str),
                   genes=q.var_names.to_numpy(dtype=str))
        info["shape"] = list(expr.shape)
        info["sum_counts"] = float(expr.sum(dtype=np.float64))
    info["wall_seconds"] = time.monotonic() - started
    info["parent_peak_rss_mib"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
    info["largest_child_peak_rss_mib"] = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
    atomic_json(out / "run.json", info)
    print(json.dumps(info), flush=True)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=["mender", "native", "cpp", "full"])
    p.add_argument("target")
    p.add_argument("tag")
    p.add_argument("--descriptor-tag", default="mender4")
    p.add_argument("--workers", type=int, default=4)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--seed", type=int, default=42)
    run(p.parse_args())
