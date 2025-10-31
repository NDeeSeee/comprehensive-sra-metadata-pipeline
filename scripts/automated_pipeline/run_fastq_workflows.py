#!/usr/bin/env python3
import argparse
import concurrent.futures as futures
import os
import shlex
import subprocess
import sys
import threading
import time
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

DEFAULT_ROOTS = [
    "/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/Controls",
    "/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/Tumors",
    "/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/Premalignant",
]

def find_cancer_dirs(roots: List[str], only_missing: bool) -> List[Path]:
    targets: List[Path] = []
    for root in roots:
        root_path = Path(root)
        if not root_path.exists():
            continue
        for sample_file in root_path.rglob("sample_list.txt"):
            cancer_dir = sample_file.parent
            status_file = cancer_dir / "sample_list.with_status.txt"
            if only_missing and status_file.exists():
                continue
            targets.append(cancer_dir)
    return sorted(set(targets), key=str)

def run_one(cancer_dir: Path, extra_args: List[str], dry_run: bool) -> Tuple[Path, int, str]:
    # Call the underlying Python workflow script directly, so we don't rely on
    # shell functions/aliases which are not visible to subprocess.
    manual_dir_env = os.environ.get("POSEIDON_MANUAL_DIR")
    if manual_dir_env:
        script_path = Path(manual_dir_env) / "sample_list_to_fastq_workflow.py"
    else:
        # Derive manual_pipeline dir relative to this file: scripts/automated_pipeline/..
        script_path = Path(__file__).resolve().parents[1] / "manual_pipeline" / "sample_list_to_fastq_workflow.py"

    cmd = ["python3", str(script_path), str(cancer_dir), *extra_args]
    if dry_run:
        return cancer_dir, 0, f"DRY-RUN: {' '.join(shlex.quote(x) for x in cmd)}"
    # Stream child output live for better observability
    try:
        print(f"[START] {cancer_dir}", flush=True)
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        # Prefix each line with directory name for clarity
        prefix = f"[{cancer_dir.name}] "
        streamed_any = False
        if proc.stdout is not None:
            for line in proc.stdout:
                streamed_any = True
                # Avoid double newlines
                sys.stdout.write(prefix + line)
                sys.stdout.flush()
        rc = proc.wait()
        # Do not duplicate output; return empty body if we streamed
        return cancer_dir, rc, ("" if streamed_any else None) or ""
    except Exception as e:
        return cancer_dir, 1, f"ERROR: {e}"

def main():
    parser = argparse.ArgumentParser(
        description="Batch-run fastq-workflow over POSEIDON subdirectories."
    )
    parser.add_argument(
        "roots",
        nargs="*",
        default=DEFAULT_ROOTS,
        help="Root directories to scan (default: Controls, Tumors, Premalignant).",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=min(8, os.cpu_count() or 4),
        help="Max parallel submissions (default: <=8).",
    )
    parser.add_argument(
        "--only-missing",
        action="store_true",
        default=True,
        help="Only run where sample_list.with_status.txt is missing (default: true).",
    )
    parser.add_argument(
        "--include-completed",
        dest="only_missing",
        action="store_false",
        help="Also run in directories that already have a status file.",
    )
    parser.add_argument(
        "--no-wait",
        action="store_true",
        default=True,
        help="Pass --no-wait to fastq-workflow (default: true).",
    )
    parser.add_argument(
        "--wait",
        dest="no_wait",
        action="store_false",
        help="Do not pass --no-wait; wait for jobs if tool supports it.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print planned commands without executing.",
    )
    parser.add_argument(
        "--status-interval",
        type=int,
        default=120,
        help="Seconds between periodic overall-status updates for the provided roots (0 disables). Default: 120.",
    )
    parser.add_argument(
        "--extra-arg",
        action="append",
        default=[],
        help="Extra argument(s) to pass through to fastq-workflow. Repeatable.",
    )
    args = parser.parse_args()

    extra_args = []
    if args.no_wait:
        extra_args.append("--no-wait")
    extra_args.extend(args.extra_arg)

    targets = find_cancer_dirs(args.roots, only_missing=args.only_missing)
    if not targets:
        print("No targets found. Check roots or flags.")
        return

    print(f"Found {len(targets)} target directories.")
    for d in targets:
        print(f"- {d}")

    # Periodic overall-status for provided roots (not global), in a background thread
    interval = int(max(0, args.status_interval or 0))
    if 0 < interval < 30:
        interval = 30  # enforce a reasonable minimum

    stop_event = threading.Event()

    def _status_script_path() -> Path:
        manual_dir_env = os.environ.get("POSEIDON_MANUAL_DIR")
        if manual_dir_env:
            return Path(manual_dir_env) / "compute_overall_status.py"
        return Path(__file__).resolve().parents[1] / "manual_pipeline" / "compute_overall_status.py"

    def _run_overall_status_once():
        script_path = _status_script_path()
        for root in args.roots:
            out_md = Path(root) / "overall_status.md"
            try:
                subprocess.run(["python3", str(script_path), str(root), str(out_md)],
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
                now = datetime.now().strftime("%H:%M:%S")
                print(f"[STATUS] updated {out_md} at {now}")
            except Exception as e:
                print(f"[STATUS] error updating for {root}: {e}")

    def _status_loop():
        if interval == 0:
            return
        # Initial run
        _run_overall_status_once()
        while not stop_event.wait(interval):
            _run_overall_status_once()

    successes, failures = 0, 0
    t_status = None
    if interval > 0:
        t_status = threading.Thread(target=_status_loop, daemon=True)
        t_status.start()

    with futures.ThreadPoolExecutor(max_workers=args.max_workers) as pool:
        jobs = {pool.submit(run_one, d, extra_args, args.dry_run): d for d in targets}
        for future in futures.as_completed(jobs):
            d, rc, out = future.result()
            prefix = "OK " if rc == 0 else "ERR"
            print(f"\n[{prefix}] {d} (rc={rc})")
            if out:
                print(out)
            if rc == 0:
                successes += 1
            else:
                failures += 1

    # Stop periodic status and do a final update
    if t_status is not None:
        stop_event.set()
        t_status.join(timeout=5)
        _run_overall_status_once()

    print(f"\nDone. Success: {successes}, Failures: {failures}")

if __name__ == "__main__":
    main()
