"""Lease the smallest available pair from a shared, resumable two-stage queue."""
import argparse
import fcntl
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time
from layers import target_layers

ROOT = Path(__file__).resolve().parent
FULL = ROOT / "full"
LIMIT = 1_200_000_000_000
HEADROOM = 50_000_000_000
MEMORY_POLL_SECONDS = 1
MAX_ATTEMPTS = 3
HOST = os.uname().nodename
JOB = os.environ["SLURM_JOB_ID"]
STEP = os.environ["SLURM_STEP_ID"]


def atomic(path, data):
    tmp = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    tmp.write_text(json.dumps(data, indent=2) + "\n")
    tmp.replace(path)


def read(path):
    return json.loads(path.read_text()) if path.exists() else {}


def node_used():
    with open("/proc/meminfo") as stream:
        info = {line.split(":")[0]: int(line.split()[1]) * 1024 for line in stream}
    return info["MemTotal"] - info["MemAvailable"]


def charged_memory():
    # At this point the supervisor and all its children are in one Slurm task cgroup.
    rel = Path("/proc/self/cgroup").read_text().split("0::", 1)[1].strip()
    cg = Path("/sys/fs/cgroup") / rel.lstrip("/")
    stats = dict(line.split() for line in (cg / "memory.stat").read_text().splitlines())
    return int((cg / "memory.current").read_text()) - int(stats.get("inactive_file", 0))


def descendants(pid):
    yield pid
    path = Path(f"/proc/{pid}/task/{pid}/children")
    try:
        children = path.read_text().split()
    except (FileNotFoundError, ProcessLookupError):
        return
    for child in children:
        yield from descendants(int(child))


def pss(pid):
    total = 0
    for member in descendants(pid):
        try:
            for line in Path(f"/proc/{member}/smaps_rollup").read_text().splitlines():
                if line.startswith("Pss:"):
                    total += int(line.split()[1]) * 1024
                    break
        except (FileNotFoundError, ProcessLookupError, OSError):
            continue
    return total


def interrupt(item, sig):
    try:
        os.killpg(item["process"].pid, sig)
    except ProcessLookupError:
        pass


def eligible(pair, stage):
    if stage == "cpu":
        key, _, _ = target_layers(pair)[0]
        return (FULL / "geometry" / pair["id"] / f"{key}.npz").is_file()
    return all((FULL / "sections" / f'{pair[k]["id"]}.h5ad').exists() for k in ("left", "right"))


def state_path(stage, ident):
    return FULL / "queue_state" / f"{stage}-{ident}.json"


def done(stage, ident):
    return (FULL / f"{stage}_done" / f"{ident}.json").exists()


def claim(pair, stage, env):
    ident = pair["id"]
    path = state_path(stage, ident)
    lock = path.with_suffix(".lock").open("a")
    try:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError:
        lock.close()
        return None
    old = read(path)
    attempts = old.get("attempt", 0)
    if done(stage, ident) or attempts >= MAX_ATTEMPTS or old.get("retry_after", 0) > time.time():
        lock.close()
        return None
    attempt = attempts + 1
    log = (FULL if FULL != ROOT / "full" else ROOT) / "logs" / f"full-{stage}-{ident}-a{attempt}.log"
    with log.open("a") as output:
        proc = subprocess.Popen([sys.executable, "-u", str(ROOT / "full_pair.py"),
                                 stage, ident, "--root", str(FULL)],
                                stdout=output, stderr=subprocess.STDOUT, env=env, start_new_session=True)
    data = {"state": "RUNNING", "stage": stage, "pair": ident, "attempt": attempt,
            "pid": proc.pid, "node": HOST, "step": f"{JOB}.{STEP}", "start_epoch": time.time(),
            "cells": pair["cells"], "targets": pair["targets"], "log": str(log)}
    atomic(path, data)
    print("CLAIM", json.dumps(data), flush=True)
    return {"process": proc, "lock": lock, "path": path, "data": data}


def retire(item, code):
    data = item["data"]
    success = done(data["stage"], data["pair"])
    data.update(exit_code=code, end_epoch=time.time())
    data["state"] = "DONE" if success else "FAILED" if data["attempt"] >= MAX_ATTEMPTS else "RETRY"
    if data["state"] == "RETRY":
        data["retry_after"] = time.time() + (120 if data.get("kill_reason") else 20)
    atomic(item["path"], data)
    item["lock"].close()
    print("RETIRE", json.dumps(data), flush=True)


def run(stage, workers):
    pairs = read(FULL / "QUEUE.json")["pairs"]
    env = os.environ.copy()
    if stage == "cpu":
        env["CUDA_VISIBLE_DEVICES"] = ""
    active = []
    if stage == "cpu":
        print("MEMORY_GUARD", LIMIT, charged_memory(), node_used(), flush=True)
    heartbeat = FULL / "queue_state" / f"worker-{stage}-{JOB}.{STEP}.json"
    last_beat = 0
    try:
        while True:
            next_poll = time.monotonic() + (MEMORY_POLL_SECONDS if stage == "cpu" else 2)
            for item in active[:]:
                code = item["process"].poll()
                if code is not None:
                    retire(item, code)
                    active.remove(item)
            used = max(node_used(), charged_memory()) if stage == "cpu" else 0
            if used > LIMIT and active and not any(i.get("killing") for i in active):
                victim = max(active, key=lambda i: pss(i["process"].pid))
                victim["data"]["kill_reason"] = "1.2 TB node memory"
                victim["data"]["node_used_bytes"] = used
                victim["data"]["task_pss_bytes"] = pss(victim["process"].pid)
                victim["killing"] = time.monotonic()
                atomic(victim["path"], victim["data"])
                print("MEMORY_INTERRUPT", json.dumps(victim["data"]), flush=True)
                interrupt(victim, signal.SIGTERM)
            for item in active:
                if item.get("killing") and time.monotonic() - item["killing"] > 4:
                    interrupt(item, signal.SIGKILL)
            if len(active) < workers and used < LIMIT - HEADROOM:
                for pair in pairs:
                    if len(active) >= workers:
                        break
                    ident = pair["id"]
                    if done(stage, ident) or not eligible(pair, stage):
                        continue
                    item = claim(pair, stage, env)
                    if item:
                        active.append(item)
            states = [{"id": p["id"], "state": "DONE" if done(stage, p["id"]) else
                       read(state_path(stage, p["id"])).get("state", "WAITING")} for p in pairs]
            if time.monotonic() - last_beat > 15:
                atomic(heartbeat, {"stage": stage, "job": JOB, "step": STEP, "node": HOST,
                                   "epoch": time.time(), "node_used_bytes": used,
                                   "active": [i["data"]["pair"] for i in active],
                                   "done": sum(s["state"] == "DONE" for s in states)})
                last_beat = time.monotonic()
            complete = all(item["state"] in {"DONE", "FAILED"} for item in states)
            if stage == "cpu" and not active and not complete:
                complete = all(s["state"] in {"DONE", "FAILED"} or
                               read(state_path("gpu", s["id"])).get("state") == "FAILED" for s in states)
            if complete and not active:
                break
            time.sleep(max(0, next_poll - time.monotonic()))
    finally:
        for item in active:
            interrupt(item, signal.SIGKILL)
            retire(item, item["process"].wait())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("stage", choices=["gpu", "cpu"])
    parser.add_argument("workers", type=int)
    parser.add_argument("--root", type=Path, default=FULL)
    args = parser.parse_args()
    FULL = args.root
    signal.signal(signal.SIGTERM, lambda *_: sys.exit(143))
    run(args.stage, args.workers)
