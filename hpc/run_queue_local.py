"""Local fallback queue runner (user policy 2026-06-12: when the HPC is
occupied/unreachable, run on the local PC keeping a 20% CPU reserve).

Runs the same hpc/manifest.txt with ceil(0.8 * ncpu) concurrent
single-threaded workers. Resume-safe: cells with an existing summary.json
are skipped, so a later HPC deployment can pick up only the remainder.

Usage:  nohup python3 hpc/run_queue_local.py > results_local_queue/queue.log 2>&1 &
"""

import os, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MANIFEST = os.path.join(ROOT, "hpc", "manifest.txt")
OUTBASE = os.path.join(ROOT, "results_hpc")   # same layout as remote
LOGDIR = os.path.join(ROOT, "results_local_queue", "logs")
WORKERS = max(1, int((os.cpu_count() or 4) * 0.8))

ENV = dict(os.environ,
           PFC_FFT_THREADS="1", OMP_NUM_THREADS="1",
           OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           VECLIB_MAXIMUM_THREADS="1")


def run_one(line):
    name, args = line.split("|", 1)
    out_dir = os.path.join(ROOT, args.split("--out ")[1].strip())
    if os.path.exists(os.path.join(out_dir, "summary.json")):
        return f"skip {name} (done)"
    t0 = time.time()
    log = os.path.join(LOGDIR, f"{name}.log")
    with open(log, "w") as lf:
        rc = subprocess.run(
            [sys.executable, os.path.join(ROOT, "hpc", "hpc_run_one.py")]
            + args.split(),
            cwd=ROOT, env=ENV, stdout=lf, stderr=subprocess.STDOUT).returncode
    return f"end {name} rc={rc} ({time.time()-t0:.0f}s)"


def main():
    os.makedirs(LOGDIR, exist_ok=True)
    os.makedirs(OUTBASE, exist_ok=True)
    lines = [l.strip() for l in open(MANIFEST) if l.strip()]
    print(f"local queue: {len(lines)} cells, {WORKERS} workers "
          f"(of {os.cpu_count()} cores, 20% reserved)", flush=True)
    t0 = time.time()
    done = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as ex:
        for msg in ex.map(run_one, lines):
            done += 1
            print(f"[{done}/{len(lines)}] {msg}", flush=True)
    print(f"queue done in {(time.time()-t0)/60:.1f} min", flush=True)


if __name__ == "__main__":
    main()
