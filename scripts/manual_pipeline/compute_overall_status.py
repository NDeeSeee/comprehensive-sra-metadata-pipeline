#!/usr/bin/env python3
import os
import sys
from collections import defaultdict

ROOT = "/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON"


def find_sample_dirs(root: str):
    sample_dirs = set()
    for dirpath, dirnames, filenames in os.walk(root):
        if "sample_list.txt" in filenames:
            sample_dirs.add(dirpath)
    return sample_dirs


def find_status_dirs(root: str):
    status_dirs = set()
    for dirpath, dirnames, filenames in os.walk(root):
        if "sample_list.with_status.txt" in filenames:
            status_dirs.add(dirpath)
    return status_dirs


def cancer_label(dir_path: str) -> str:
    # Use first two path components under ROOT as the cancer type label
    rel = os.path.relpath(dir_path, ROOT)
    parts = rel.split(os.sep)
    return "/".join(parts[:2]) if len(parts) >= 2 else parts[0]


def tally_status(status: str, acc: dict):
    if status == "DBGaP_REQUIRED":
        acc["dbgap_required"] += 1
    elif status == "BAM_DONE":
        acc["bam_done"] += 1
    elif status == "FASTQ_DONE":
        acc["fastq_done"] += 1
    elif status.startswith("PENDING (FASTQ -> BAM"):
        acc["pending_fastq_to_bam"] += 1
    elif status.startswith("PENDING (SRA -> FASTQ"):
        acc["pending_sra_to_fastq"] += 1
    elif status.startswith("PENDING (SRA downloading"):
        acc["pending_sra_downloading"] += 1
    elif status == "PENDING":
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
    if done_fastq == exp and exp > 0 and done_bam < exp:
        return "OPTIONAL_BAM_ALIGNMENT"
    if exp > (blocked + done_fastq + done_bam):
        return "START_SRA_DOWNLOADING"
    if blocked > 0 and (done_bam + done_fastq + blocked) == exp:
        return "REQUEST_ACCESS"
    return "CHECK_MANUALLY"


def main():
    # Usage: compute_overall_status.py [root_dir] [out_md]
    root = sys.argv[1] if len(sys.argv) > 1 else ROOT
    out_path = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "Master_Project", "overall_status.md")

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
        per_dir_stats[dir_path] = (cancer_label(dir_path), action_required, acc)

    # Build Markdown table content
    header = "| " + " | ".join(essential_cols) + " |"
    sep_cells = ["---"] + ["---:"] * (len(essential_cols) - 1)
    separator = "| " + " | ".join(sep_cells) + " |"
    lines = [header, separator]

    for dir_path in sorted(per_dir_stats.keys(), key=lambda p: cancer_label(p)):
        label, action_required, acc = per_dir_stats[dir_path]
        if action_required.startswith("RUN "):
            # Replace numeric columns with None for RUN actions to avoid misleading zeros
            cells = [label, action_required] + ["None"] * 8
        else:
            cells = [
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
        lines.append("| " + " | ".join(cells) + " |")

    # Ensure output directory exists and write the table (.md)
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")

    # Also print path to the saved table
    print(out_path)


if __name__ == "__main__":
    main()
