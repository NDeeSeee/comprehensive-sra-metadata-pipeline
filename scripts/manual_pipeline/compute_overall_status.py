#!/usr/bin/env python3
import os
import sys
import csv
from collections import defaultdict

# Default project root; allow override via CLI arg or POSEIDON_DIR env var
DEFAULT_ROOT = "/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON"


def _prune_heavy_dirs(dirnames: list):
    """Remove directories we never need to descend into for status discovery."""
    # Keep this conservative; only skip well-known heavy dirs
    skip = {
        "Master_Project",  # contains outputs like overall_status.md
        "star_output",     # aggregated STAR outputs, no sample_list files
        "logs",            # log directories can be huge
        "bams",            # BAM outputs, not holding sample_list files
        ".git",
        "__pycache__",
    }
    # Modify in-place
    dirnames[:] = [d for d in dirnames if d not in skip]


def find_sample_dirs(root: str):
    sample_dirs = set()
    for dirpath, dirnames, filenames in os.walk(root, topdown=True):
        _prune_heavy_dirs(dirnames)
        if "sample_list.txt" in filenames:
            sample_dirs.add(dirpath)
    return sample_dirs


def find_status_dirs(root: str):
    status_dirs = set()
    for dirpath, dirnames, filenames in os.walk(root, topdown=True):
        _prune_heavy_dirs(dirnames)
        if "sample_list.with_status.txt" in filenames:
            status_dirs.add(dirpath)
    return status_dirs


def cancer_label(dir_path: str, root: str) -> str:
    """Return a concise label relative to the effective root.

    Prefer the first two path components to keep the table width reasonable,
    but work correctly for any root passed at runtime.
    """
    rel = os.path.relpath(dir_path, root)
    if rel == ".":  # if the directory is the root itself, show its name
        return os.path.basename(os.path.abspath(root))
    parts = rel.split(os.sep)
    return "/".join(parts[:2]) if len(parts) >= 2 else parts[0]


def tally_status(status: str, acc: dict):
    # Normalize to new status names from sample_list_to_fastq_workflow
    if status == "DBGaP_REQUIRED":
        acc["dbgap_required"] += 1
    elif status == "BAM_DONE":
        acc["bam_done"] += 1
    elif status == "FASTQ_DONE":
        acc["fastq_done"] += 1
    elif status == "ALIGN_IN_PROGRESS":
        acc["pending_fastq_to_bam"] += 1
    elif status.startswith("CONVERTING"):
        acc["pending_sra_to_fastq"] += 1
    elif status == "NEEDS_CONVERSION":
        acc["pending_sra_to_fastq"] += 1
    elif status == "NEEDS_PREFETCH":
        acc["pending_sra_downloading"] += 1
    elif status in ("PENDING", "UNKNOWN"):
        acc["pending"] += 1
    else:
        acc["pending"] += 1  # fallback


essential_cols = [
    "cancer_type",
    "action_required",
    "expected",
    "dbgap_required",
    "pending_sra_downloading",
    "pending_sra_to_fastq",
    "pending_fastq_to_bam",
    "pending",
    "fastq_done",
    "bam_done",
]


def compute_action_required(has_sample_list: bool, has_status: bool, acc: dict) -> str:
    # 3a
    if not has_sample_list:
        return "RUN generate_geo_sample_lists.py (missing sample_list.txt)"
    # 3b
    if has_sample_list and not has_status:
        return "RUN sample_list_to_fastq_workflow.py (missing sample_list.with_status.txt)"
    # 3c with tallies
    exp = acc.get("expected", 0)
    done_fastq = acc.get("fastq_done", 0)
    done_bam = acc.get("bam_done", 0)
    blocked = acc.get("dbgap_required", 0)
    any_pending = (
        acc.get("pending_sra_downloading", 0)
        + acc.get("pending_sra_to_fastq", 0)
        + acc.get("pending_fastq_to_bam", 0)
        + acc.get("pending", 0)
    ) > 0

    # Priority order
    if any_pending:
        return "WAIT"
    if done_bam == exp and exp > 0:
        return "NONE"
    # optional alignment when all samples are at least FASTQ_DONE or BAM_DONE
    if (done_bam + done_fastq) == exp and exp > 0 and done_bam < exp:
        return "OPTIONAL_BAM_ALIGNMENT"
    if exp > (blocked + done_fastq + done_bam):
        return "START_SRA_DOWNLOADING"
    # Do not emit REQUEST_ACCESS here; we already show dbgap_required in a dedicated column
    if blocked > 0 and (done_bam + done_fastq + blocked) == exp:
        return "NONE"
    return "CHECK_MANUALLY"


def main():
    # Usage: compute_overall_status.py [root_dir] [out_md]
    if len(sys.argv) > 1:
        root = os.path.abspath(sys.argv[1])
    else:
        root = os.path.abspath(os.environ.get("POSEIDON_DIR", DEFAULT_ROOT))

    out_path = (
        sys.argv[2]
        if len(sys.argv) > 2
        else os.path.join(root, "Master_Project", "overall_status.md")
    )

    sample_dirs = find_sample_dirs(root)
    status_dirs = find_status_dirs(root)
    all_dirs = sorted(sample_dirs | status_dirs)

    per_dir_stats = {}

    for dir_path in all_dirs:
        has_sample_list = os.path.isfile(os.path.join(dir_path, "sample_list.txt"))
        status_path = os.path.join(dir_path, "sample_list.with_status.txt")
        has_status = os.path.isfile(status_path)

        acc = defaultdict(int)
        if has_status:
            # Parse status file; support 5-col (new) and 4-col (legacy) formats
            try:
                with open(status_path, "r", encoding="utf-8", errors="ignore") as f:
                    for line in f:
                        s = line.strip()
                        if not s:
                            continue
                        parts = s.split("\t")
                        status = None
                        if len(parts) == 5:
                            # SAMPLE_ID, ACTION_REQUIRED, R1_LIST, R2_LIST, STATUS
                            status = parts[4]
                        elif len(parts) == 4:
                            # Legacy: SAMPLE_ID, R1_LIST, R2_LIST, STATUS
                            status = parts[3]
                        else:
                            continue
                        acc["expected"] += 1
                        tally_status(status, acc)
            except FileNotFoundError:
                pass
        # Compute directory-level action_required
        action_required = compute_action_required(has_sample_list, has_status, acc)
        per_dir_stats[dir_path] = (cancer_label(dir_path, root), action_required, acc)

    # Prepare rows in a common structure
    rows = []
    for dir_path in sorted(per_dir_stats.keys(), key=lambda p: cancer_label(p, root)):
        label, action_required, acc = per_dir_stats[dir_path]
        if action_required.startswith("RUN "):
            row = [label, action_required] + [""] * 8
        else:
            row = [
                label,
                action_required,
                str(acc.get("expected", 0)),
                str(acc.get("dbgap_required", 0)),
                str(acc.get("pending_sra_downloading", 0)),
                str(acc.get("pending_sra_to_fastq", 0)),
                str(acc.get("pending_fastq_to_bam", 0)),
                str(acc.get("pending", 0)),
                str(acc.get("fastq_done", 0)),
                str(acc.get("bam_done", 0)),
            ]
        rows.append(row)

    # Write CSV if requested, otherwise write Markdown for backward compatibility
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    if out_path.lower().endswith(".csv"):
        with open(out_path, "w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(essential_cols)
            writer.writerows(rows)
    else:
        header = "| " + " | ".join(essential_cols) + " |"
        sep_cells = [":---" if i < 2 else "---:" for i in range(len(essential_cols))]
        separator = "| " + " | ".join(sep_cells) + " |"
        lines = [header, separator]
        for cells in rows:
            lines.append("| " + " | ".join(cells) + " |")
        with open(out_path, "w", encoding="utf-8") as f:
            f.write("\n".join(lines) + "\n")

    print(out_path)


if __name__ == "__main__":
    main()
