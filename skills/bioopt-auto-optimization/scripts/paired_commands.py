"""Run resource-matched command pairs and retain a validation receipt for each."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import signal
import statistics
import subprocess
import time
import threading

THREAD_VARIABLES = (
    "OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS",
)
COMMAND_KEYS = ("baseline", "candidate", "validate")


def file_digest(path):
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for chunk in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def expand(argv, values):
    for value in argv:
        for name, replacement in values.items():
            value = value.replace("{" + name + "}", str(replacement))
        yield value


def run(argv, folder, label, env, timeout):
    start = time.monotonic()
    timed_out = threading.Event()
    with (folder / f"{label}.log").open("wb") as log:
        try:
            proc = subprocess.Popen(
                argv, stdout=log, stderr=subprocess.STDOUT, env=env, start_new_session=True
            )
        except OSError as exc:
            log.write(f"{exc}\n".encode())
            return {"status": "spawn_failed", "elapsed_seconds": None, "peak_rss_kib": None}

        def kill_process_group():
            try:
                os.killpg(proc.pid, signal.SIGKILL)
                timed_out.set()
            except ProcessLookupError:
                pass

        timer = threading.Timer(timeout, kill_process_group)
        timer.start()
        try:
            _, status, usage = os.wait4(proc.pid, 0)
        finally:
            timer.cancel()
            timer.join()
        proc.returncode = os.waitstatus_to_exitcode(status)
    return {
        "elapsed_seconds": time.monotonic() - start,
        "peak_rss_kib": usage.ru_maxrss,
        "status": "timeout" if timed_out.is_set() else proc.returncode,
    }


def save(report, path):
    sibling = path.with_suffix(".tmp")
    sibling.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    sibling.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("manifest", type=Path)
    parser.add_argument("--out", type=Path, required=True, help="A new evidence directory")
    args = parser.parse_args()
    spec = json.loads(args.manifest.read_text())
    try:
        for key in COMMAND_KEYS:
            argv = spec[key]
            if not isinstance(argv, list) or not argv or not all(isinstance(a, str) for a in argv):
                raise ValueError(f"{key} must be a nonempty argv list of strings")
        case_id = spec["case_id"]
        scope = spec["claim_scope"]
        contract = spec["equivalence_contract"]
        repeats = int(spec.get("repeats", 5))
        warmups = int(spec.get("warmups", 1))
        timeout = int(spec.get("timeout_seconds", 300))
        cpu_ids = sorted(os.sched_getaffinity(0))
        threads = int(spec.get("threads", len(cpu_ids)))
        if not scope or not contract or repeats < 1 or warmups < 0 or timeout < 1:
            raise ValueError("supply a scope and contract, with positive repeats and timeout")
        if not 1 <= threads <= len(cpu_ids):
            raise ValueError(f"threads must fit the {len(cpu_ids)} CPU affinity slots")
    except (ValueError, TypeError, KeyError) as exc:
        parser.error(str(exc))

    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    env_changes = dict.fromkeys(THREAD_VARIABLES, str(threads))
    env = {**os.environ, **env_changes}
    receipt = {
        "schema_version": "1.0",
        "case_id": case_id,
        "claim_scope": scope,
        "contract": contract,
        "input_sha256": {key: file_digest(Path(key)) for key in spec.get("inputs", [])},
        "environment": {
            "host": os.uname().nodename,
            "cpu_ids": cpu_ids,
            "threads": env_changes,
            "slurm_job": os.environ.get("SLURM_JOB_ID"),
            "slurm_step": os.environ.get("SLURM_STEP_ID"),
            "cuda_visible_devices": os.environ.get("CUDA_VISIBLE_DEVICES"),
        },
        "commands": {key: spec[key] for key in COMMAND_KEYS},
        "warmups": warmups,
        "repeats": repeats,
        "pairs": [],
        "status": "incomplete",
    }

    for index in range(warmups + repeats):
        folder = out / f"pair_{index:03d}"
        folder.mkdir()
        paths = {key: folder / key for key in COMMAND_KEYS[:2]}
        for path in paths.values():
            path.mkdir()
        order = list(COMMAND_KEYS[:2])
        if index % 2:
            order.reverse()
        pair = {"warmup": index < warmups, "order": order, "runs": {}}
        receipt["pairs"].append(pair)
        for label in [*order, "validate"]:
            values = {**paths, "pair": index, "output": paths.get(label, folder)}
            cmd = list(expand(spec[label], values))
            scratch = folder / f"{label}_tmp"
            scratch.mkdir()
            pair["runs"][label] = run(cmd, folder, label, {**env, "TMPDIR": str(scratch)}, timeout)
            if pair["runs"][label]["status"] != 0:
                receipt["status"] = "validation_failed" if label == "validate" else "run_failed"
                save(receipt, out / "receipt.json")
                return 1
        save(receipt, out / "receipt.json")

    measured = [item["runs"] for item in receipt["pairs"] if not item["warmup"]]
    paired_ratios = [
        row["baseline"]["elapsed_seconds"] / row["candidate"]["elapsed_seconds"] for row in measured
    ]
    receipt["summary"] = {
        label: {
            "median_seconds": statistics.median(row[label]["elapsed_seconds"] for row in measured),
            "range_seconds": [min(row[label]["elapsed_seconds"] for row in measured),
                              max(row[label]["elapsed_seconds"] for row in measured)],
            "max_observed_rss_kib": max(
                (row[label]["peak_rss_kib"] for row in measured if row[label]["peak_rss_kib"] is not None),
                default=None,
            ),
        } for label in COMMAND_KEYS[:2]
    }
    receipt["summary"]["median_paired_speedup"] = statistics.median(paired_ratios)
    receipt["status"] = "pass"
    save(receipt, out / "receipt.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
