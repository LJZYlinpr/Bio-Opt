import argparse, hashlib, json, os, shutil, subprocess, sys, time
from pathlib import Path

parser = argparse.ArgumentParser(description="Repeat UniST interpolation with 2, 4, and 8 threads")
parser.add_argument("--source", type=Path, required=True, help="Pinned UniST checkout")
parser.add_argument("--out", type=Path, required=True, help="New working/output directory")
args = parser.parse_args()
SRC = args.source.resolve()
W = args.out.resolve()
W.mkdir(parents=True, exist_ok=False)
PY = sys.executable
OUT = W


def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def run(threads, tag):
    wd = W / f"repeat_{tag}"
    (wd / "input5").mkdir(parents=True)
    shutil.copy(SRC / "example_interpolation/slice1.tif", wd / "input5/slice1.tif")
    shutil.copy(SRC / "example_interpolation/slice5.tif", wd / "input5/slice5.tif")
    env = os.environ.copy()
    env.update(PYTHONPATH=str(SRC), MPLCONFIGDIR=str(W / "mplcache"), TF_CPP_MIN_LOG_LEVEL="1",
               CUDA_VISIBLE_DEVICES="", OMP_NUM_THREADS=str(threads),
               TF_NUM_INTRAOP_THREADS=str(threads), TF_NUM_INTEROP_THREADS=str(threads))
    steps = [
        ("convert", [PY, "-m", "interpolation.convert_images", "--input_dir", "input5", "--output_dir", "input5_3ch"]),
        ("interp", [PY, "-m", "interpolation.main", "--mode", "auto", "--tile_size", "1024", "1024", "--pth", "input5_3ch/"]),
        ("merge", [PY, "-m", "interpolation.merge_volume", "--base_dir", "input5", "--interp_dir", "input5_3ch",
                   "--out_binary_dir", "binary_slices", "--out_stack_path", "volume.tif", "--threshold", "127"]),
    ]
    times = {}
    t_all = time.time()
    for name, cmd in steps:
        t0 = time.time()
        p = subprocess.run(cmd, cwd=wd, env=env, capture_output=True, text=True)
        times[name] = round(time.time() - t0, 3)
        if p.returncode != 0:
            raise SystemExit(f"{tag}/{name} rc={p.returncode}\n{p.stderr[-500:]}")
    return {"threads": threads, "tag": tag, "total_s": round(time.time() - t_all, 3),
            "volume_sha": sha(wd / "volume.tif"), "steps": times}


seq = [2, 4, 8, 2, 4, 8]
runs = []
for i, t in enumerate(seq):
    r = run(t, f"{t}_r{i//3}")
    runs.append(r)
    print(json.dumps(r), flush=True)

shas = {r["volume_sha"] for r in runs}
summary = {"sequence": seq, "all_identical": len(shas) == 1, "unique_sha": sorted(shas),
           "runs": runs,
           "mean_by_threads": {t: round(sum(r["total_s"] for r in runs if r["threads"] == t) / sum(1 for r in runs if r["threads"] == t), 3) for t in sorted(set(seq))}}
OUT.mkdir(parents=True, exist_ok=True)
(OUT / "unist_threads_repeat_RESULT.json").write_text(json.dumps(summary, indent=2) + "\n")
print(json.dumps(summary["mean_by_threads"]), flush=True)
print("ALL_IDENTICAL", summary["all_identical"], flush=True)

if len(shas) != 1:
    raise SystemExit("UniST output equivalence failed")
