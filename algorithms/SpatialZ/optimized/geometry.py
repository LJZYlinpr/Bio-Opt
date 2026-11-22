"""Run the unchanged SpatialZ geometry to the verified CPU boundary."""
import argparse
import importlib.util
import json
from pathlib import Path
import os
import time

import numpy as np
import torch

ROOT = Path(__file__).resolve().parent
FIG3 = ROOT.parent / "fig3_three_window_reconstruction_20261120"
NATIVE = ROOT.parent / "original/SpatialZ.py"
CONFIG = {
    "lr": 1e5, "nb_iter_max": 1000, "num_projections": 200, "seed": 42,
    "n_mag": 1.0, "k_sam": 50, "n_cell": None, "syn_mode": "default", "Beta": 100,
}


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class GeometryReady(Exception):
    pass


def write_geometry(native, loader, frames, plane, out):
    started = time.monotonic()
    out.mkdir(parents=True, exist_ok=True)
    left_id, right_id = plane["left_visible"], plane["right_visible"]
    a, b = frames[left_id].copy(), frames[right_id].copy()
    for frame in (a, b):
        frame.X = loader.dense(frame.X)
        frame.obsm["spatial"] = np.asarray(frame.obsm["spatial"])[:, :2]
    captured = []

    def capture(a, b, query, label, **kwargs):
        assert label == "cell_class" and kwargs == {
            "scale": 6, "radius": 15, "micro_env_key": "mender"
        }
        captured.extend((a, b, query))
        raise GeometryReady

    native.perform_microenv_measure = capture
    cuda_start = time.monotonic()
    try:
        native.Generate_spatialz(
            a, b, adata1_id=left_id, adata2_id=right_id,
            alpha=1 - plane["relative_depth"], device="cuda:0",
            cell_type_key="cell_class", verbose=True, **CONFIG,
        )
    except GeometryReady:
        pass
    cuda_wall = time.monotonic() - cuda_start
    assert len(captured) == 3
    for name, frame in zip(("left", "right", "query"), captured):
        frame.write_h5ad(out / f"{name}.h5ad")
    info = {
        "plane": plane, "config": CONFIG, "wall_to_cpu_boundary_seconds": cuda_wall,
        "wall_with_checkpoint_seconds": time.monotonic() - started,
        "rows": [f.n_obs for f in captured], "genes": captured[-1].n_vars,
        "step": os.environ.get("SLURM_STEP_ID"),
        "job": os.environ.get("SLURM_JOB_ID"),
        "gpu": str(torch.cuda.get_device_properties(0).uuid),
    }
    stage = out / "geometry.tmp"
    stage.write_text(json.dumps(info, indent=2) + "\n")
    stage.replace(out / "geometry.json")
    print(json.dumps(info), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("targets", nargs="+")
    args = parser.parse_args()
    torch.set_num_threads(1)
    assert torch.cuda.is_available()
    assert torch.ones(2, device="cuda").sum().item() == 2
    native = load_module("spatialz_native", NATIVE)
    loader = load_module("fig3_loader", FIG3 / "adapters/joint_loader.py")
    frames, _, planes, _, _ = loader.load_case("ABCA__native__FIVE20", "SpatialZ")
    lookup = {plane["section"]: plane for plane in planes}
    for target in args.targets:
        write_geometry(native, loader, frames, lookup[target], ROOT / "geometry" / target)


if __name__ == "__main__":
    main()
