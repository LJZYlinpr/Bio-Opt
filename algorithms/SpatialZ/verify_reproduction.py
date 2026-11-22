"""Check equivalence of new SpatialZ reproduction outputs."""
import argparse
import json
from pathlib import Path

import numpy as np


def require(condition, message):
    if not condition:
        raise SystemExit(f"FAIL: {message}")


def check_expression(a, b, shape, label):
    require(a.shape == b.shape == shape, f"{label}: unexpected expression shape")
    require(np.isfinite(a).all() and np.isfinite(b).all(), f"{label}: nonfinite expression")
    require(np.array_equal(a, b), f"{label}: expression values differ")


def verify_abca(directory):
    info = json.loads((directory / "comparison.json").read_text())
    require(info["exact"] is True and info["different_entries"] == 0
            and info["max_abs_error"] == 0, "ABCA 40 µm: report shows a difference")
    a = np.load(directory / "original.npy", allow_pickle=False)
    b = np.load(directory / "optimized.npy", allow_pickle=False)
    check_expression(a, b, (8756, 1122), "ABCA 40 µm")
    print(f"PASS: ABCA 40 µm, {a.size:,} expression values match.")


def verify_merfish(directory):
    info = json.loads((directory / "merfish_RESULT.json").read_text())
    expected = [(1, 2, .75, 5505), (2, 3, .5, 5741), (4, 5, .25, 5608)]
    require(len(info["cases"]) == len(expected), "MERFISH: expected exactly three cases")
    total = 0
    for i, (case, (left, right, alpha, cells)) in enumerate(zip(info["cases"], expected), 1):
        require((case["left"], case["right"], case["alpha"], case["cells"], case["genes"])
                == (left, right, alpha, cells, 155), f"MERFISH case {i}: unexpected input or output size")
        require(all(case[k] is True for k in ("xy_exact", "labels_exact", "expression_exact")) and
                case["mender_exact"] == [True, True, True] and
                case["expression_max_abs_error"] == 0 and case["expression_different_entries"] == 0,
                f"MERFISH case {i}: report shows a difference")
        with np.load(directory / f"merfish_{i}_comparison.npz", allow_pickle=False) as data:
            check_expression(data["original"], data["optimized"], (cells, 155), f"MERFISH case {i}")
            require(data["original_xy"].shape == (cells, 2) and
                    np.array_equal(data["original_xy"], data["optimized_xy"]),
                    f"MERFISH case {i}: XY differs")
            require(data["original_labels"].shape == (cells,) and
                    np.array_equal(data["original_labels"], data["optimized_labels"]),
                    f"MERFISH case {i}: cell type differs")
        total += cells * 155
    print(f"PASS: MERFISH, three cases and {total:,} expression values match.")


p = argparse.ArgumentParser(description=__doc__)
p.add_argument("--abca", type=Path, help="Output directory from reproduce_abca.py")
p.add_argument("--official", type=Path, help="Output directory from reproduce_official.py")
args = p.parse_args()
if not (args.abca or args.official):
    p.error("provide --abca, --official, or both")
if args.abca:
    verify_abca(args.abca)
if args.official:
    verify_merfish(args.official)
