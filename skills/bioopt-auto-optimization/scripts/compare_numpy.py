#!/usr/bin/env python3
"""Compare persisted NPY/NPZ fields under an explicit exact or floating contract."""
import argparse
import json
import math
from pathlib import Path

import numpy as np


def load_fields(path):
    if path.suffix == ".npy":
        return {"array": np.load(path, allow_pickle=False)}
    if path.suffix == ".npz":
        with np.load(path, allow_pickle=False) as archive:
            if len(set(archive.files)) != len(archive.files):
                raise ValueError(f"Duplicate archive fields in {path}")
            return {key: archive[key] for key in archive.files}
    raise ValueError(f"Expected .npy or .npz, got {path}")


def compare(a, b, tolerance, atol, rtol):
    if a.shape != b.shape or a.dtype != b.dtype or a.dtype.kind not in "biufcSU":
        return {"pass": False, "reason": "shape, dtype, or supported primitive type differs"}
    if tolerance and a.dtype.kind != "f":
        return {"pass": False, "reason": "tolerance requires a real floating field"}
    if a.dtype.kind in "fc" and (not np.isfinite(a).all() or not np.isfinite(b).all()):
        return {"pass": False, "reason": "non-finite values"}
    if not tolerance:
        ref_bytes = a.ravel(order="C").view(np.uint8)
        got_bytes = b.ravel(order="C").view(np.uint8)
        return {"pass": bool(np.array_equal(ref_bytes, got_bytes)), "criterion": "bitwise"}
    with np.errstate(over="ignore", invalid="ignore"):
        errors = np.abs(b.astype(np.longdouble) - a.astype(np.longdouble))
        threshold = atol + rtol * np.abs(a.astype(np.longdouble))
        exceed = int(np.count_nonzero(~np.isfinite(errors) | ~np.isfinite(threshold) | (errors > threshold)))
    return {"pass": exceed == 0, "criterion": "tolerance", "exceeding_entries": exceed,
            "max_absolute_error": str(errors.max(initial=0))}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("baseline", type=Path)
    parser.add_argument("candidate", type=Path)
    parser.add_argument("--float-key", action="append", default=[], help="Named field allowed to differ numerically; use 'array' for NPY")
    parser.add_argument("--atol", type=float, default=0.0)
    parser.add_argument("--rtol", type=float, default=0.0)
    args = parser.parse_args()
    if any(not math.isfinite(x) or x < 0 for x in (args.atol, args.rtol)):
        parser.error("tolerances must be finite and nonnegative")
    try:
        a, b = load_fields(args.baseline), load_fields(args.candidate)
        if a.keys() != b.keys() or not a:
            raise ValueError("Archive fields differ or are empty")
        if set(args.float_key) - a.keys():
            raise ValueError("Unknown tolerance field")
        fields = {key: compare(a[key], b[key], key in args.float_key, args.atol, args.rtol)
                  for key in sorted(a)}
        passed = all(field["pass"] for field in fields.values())
    except (ValueError, TypeError, OSError, OverflowError) as exc:
        print(json.dumps({"pass": False, "error": str(exc)}))
        return 2
    print(json.dumps({"pass": passed, "fields": fields, "atol": args.atol, "rtol": args.rtol}))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
