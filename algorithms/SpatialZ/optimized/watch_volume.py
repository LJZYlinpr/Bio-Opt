"""Publish progress and, after every pair succeeds, the native-volume layer index."""
import argparse
import json
import math
from pathlib import Path
from datetime import datetime
import time
from zoneinfo import ZoneInfo
from layers import target_layers

ROOT = Path(__file__).resolve().parent / "full"


def write(name, data):
    path = ROOT / name
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(path)


def read(path):
    try:
        return json.loads(path.read_text())
    except FileNotFoundError:
        return {}


def main():
    q = read(ROOT / "QUEUE.json")
    count = sum(len(p["targets"]) for p in q["pairs"])
    pairs = q["pairs"]
    while True:
        state = {"updated_cst": datetime.now(ZoneInfo("Asia/Shanghai")).isoformat(),
                 "state": "RUNNING", "z_step_um": q.get("z_step_um", 20),
                 "pairs": len(pairs), "synthetic_planes": count, "regular_planes": len(q["regular_grid_um"])}
        for stage, folder, ext in (("gpu", "geometry", "npz"), ("cpu", "predictions", "h5ad")):
            complete = {p.stem for p in (ROOT / f"{stage}_done").glob("*.json")}
            ready = sum((ROOT / folder / p["id"] / f"{key}.{ext}").is_file()
                        for p in pairs for key, _, _ in target_layers(p))
            by_status = {"RUNNING": [], "FAILED": [], "RETRY": []}
            for f in (ROOT / "queue_state").glob(f"{stage}-*.json"):
                if f.stem.removeprefix(f"{stage}-") not in complete:
                    status = read(f)
                    value = status.get("state")
                    if value in by_status:
                        by_status[value].append(status)
            state[stage] = {"pairs_done": len(complete), "planes_done": ready,
                            "active": by_status["RUNNING"], "failed": by_status["FAILED"],
                            "retry": by_status["RETRY"]}
        state["workers"] = [w for p in (ROOT / "queue_state").glob("worker-*.json")
                            if (w := read(p)) and time.time() - w["epoch"] < 45]
        write("LIVE_STATUS.json", state)
        if state["cpu"]["pairs_done"] == len(pairs):
            assert state["cpu"]["planes_done"] == count
            break
        time.sleep(20)

    import h5py
    regular = q["regular_grid_um"]
    row_info = []
    for level in q["levels"]:
        path = ROOT / "sections" / f'{level["id"]}.h5ad'
        with h5py.File(path) as h:
            assert tuple(h["X"].attrs["shape"]) == (level["cells"], q["gene_count"])
            assert h["var/_index"].asstr()[:].tolist() == q["gene_intersection"]
        row_info.append({"kind": "measured", "path": str(path), "source_sections": level["sections"],
                         "z_um": level["z_um"], "cells": level["cells"],
                         "on_regular_grid": level["id"] in q["observed_on_grid"].values()})
    for pair in pairs:
        for key, z, _ in target_layers(pair):
            path = ROOT / "predictions" / pair["id"] / f"{key}.h5ad"
            with h5py.File(path) as h:
                rows, genes = (int(n) for n in h["X"].attrs["shape"])
                assert genes == q["gene_count"] and rows == len(h["obsm/spatial"]) > 0
                assert h["var/_index"].asstr()[:].tolist() == q["gene_intersection"]
                assert h["uns/z_um"][()] == z
            row_info.append({"kind": "interpolated", "path": str(path), "source_pair": pair["id"],
                             "z_um": z, "cells": rows, "on_regular_grid": True})
    actual = sorted(r["z_um"] for r in row_info if r["on_regular_grid"])
    assert len(actual) == len(regular)
    assert all(math.isclose(x, y, rel_tol=0, abs_tol=1e-6) for x, y in zip(actual, regular))
    write("VOLUME_COMPLETE.json", {"layers": sorted(row_info, key=lambda r: r["z_um"]),
                                  "expression_scale": "raw_counts", "coordinate_units": "um",
                                  "gene_ids": q["gene_intersection"], "run_contract": str(ROOT / "QUEUE.json")})
    state.update(state="COMPLETE", updated_cst=datetime.now(ZoneInfo("Asia/Shanghai")).isoformat())
    write("LIVE_STATUS.json", state)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=ROOT)
    ROOT = parser.parse_args().root
    main()
