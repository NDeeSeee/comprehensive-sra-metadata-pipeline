#!/usr/bin/env python3

# Generate GEO sample lists (previously GEO_sampleSetup_enhanced_VP.py)
# Creates POSEIDON-compatible sample_list.txt files from Excel/CSV metadata

import os
import sys
import argparse
from collections import defaultdict
from typing import List, Optional, Dict, Tuple

import pandas as pd
import subprocess
import re
from urllib import request as urlrequest
from urllib import error as urlerror

# Hardcoded project directory - script must be run from this location
POSEIDON_DIR = "/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON"

def check_working_directory():
    """Check if script is running from the correct POSEIDON directory."""
    current_dir = os.getcwd()
    if current_dir != POSEIDON_DIR:
        print(f"ERROR: You are not in the correct directory!")
        print(f"Current directory: {current_dir}")
        print(f"Required directory: {POSEIDON_DIR}")
        print("Please go to the POSEIDON directory and run the script from there.")
        sys.exit(1)


def infer_cancer_type_from_filename(file_path: str) -> str:
    """Infer cancer type from the input filename by stripping common words.

    Example: "Pancreas metadata.xlsx" -> "Pancreas".
    Falls back to base stem if no better inference is possible.
    """
    base = os.path.splitext(os.path.basename(file_path))[0]
    # Normalize separators and split
    tokens = base.replace("-", " ").replace("_", " ").split()
    filtered_tokens = [t for t in tokens if t.lower() not in {"metadata", "meta"}]
    if not filtered_tokens:
        return base
    # Join with spaces for readability in directory names
    return " ".join(filtered_tokens).strip()


def _resolve_column_name(df: pd.DataFrame, candidates: List[str]) -> Optional[str]:
    """Return the actual column name matching any candidate (case/space-insensitive)."""
    normalized_to_actual = {str(c).strip().lower(): c for c in df.columns}
    for candidate in candidates:
        actual = normalized_to_actual.get(candidate.strip().lower())
        if actual is not None:
            return actual
    return None


def process_dataframe(df: pd.DataFrame, sheet_name: str, cancer_type: str) -> None:
    """Build BioSample -> FASTQ mapping and write sample_list.txt under
    POSEIDON/<sheet_name>/<cancer_type>/sample_list.txt.
    Requires columns akin to BioSample and Run (case-insensitive).
    """
    biosample_column = _resolve_column_name(df, ["BioSample"])
    run_column = _resolve_column_name(df, ["Run"])

    if biosample_column is None or run_column is None:
        print(
            f"Missing required columns in sheet '{sheet_name}'. "
            f"Found: {list(df.columns)} | Need: BioSample, Run"
        )
        return

    mapping = defaultdict(lambda: {"srr_ids": [], "fastq_read1": [], "fastq_read2": []})
    layout_cache: Dict[str, str] = {}

    for _, row in df.iterrows():
        biosample_raw = row.get(biosample_column)
        run_raw = row.get(run_column)

        if pd.isna(biosample_raw) or pd.isna(run_raw):
            continue

        biosample_id = str(biosample_raw).strip()
        srr_id = str(run_raw).strip()
        if not biosample_id or not srr_id:
            continue

        # Determine library layout (PAIRED or SINGLE)
        layout = _resolve_layout_for_run(srr_id, layout_cache)
        fastq_r1 = f"{srr_id}_1.fastq.gz"
        # Avoid duplicate SRR entries per BioSample
        entry = mapping[biosample_id]
        if srr_id not in entry["srr_ids"]:
            entry["srr_ids"].append(srr_id)
            entry["fastq_read1"].append(fastq_r1)
            if layout == "PAIRED":
                entry["fastq_read2"].append(f"{srr_id}_2.fastq.gz")

    # Prepare output directory and write results
    output_dir = os.path.join(POSEIDON_DIR, sheet_name, cancer_type)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "sample_list.txt")

    def _stable_unique(seq):
        seen = set()
        out = []
        for x in seq:
            if x not in seen:
                seen.add(x)
                out.append(x)
        return out

    with open(output_file, "w") as fh:
        for biosample_id, info in mapping.items():
            r1_list = _stable_unique(info["fastq_read1"])
            r2_list = _stable_unique(info["fastq_read2"])
            fastq_read1_joined = ",".join(r1_list)
            fastq_read2_joined = ",".join(r2_list) if r2_list else "NA"
            fh.write(f"{biosample_id}\t{fastq_read1_joined}\t{fastq_read2_joined}\n")

    # Optional console output for traceability
    for biosample_id, info in mapping.items():
        print(
            f"BioSample: {biosample_id} -> Associated SRR IDs: {', '.join(info['srr_ids'])}"
        )
    print(
        f"Wrote {len(mapping)} BioSample rows to: {output_file}"
    )


# CLI
parser = argparse.ArgumentParser(
    description=(
        "Generate POSEIDON sample_list.txt from metadata or SRR list. "
        "Accepts Excel (.xlsx/.xls), CSV (.csv), or a plain text file (.txt) "
        "with one SRR/ERR accession per line."
    )
)
parser.add_argument(
    "metadata_file",
    help=(
        "Path to input: Excel/CSV metadata (expects BioSample & Run columns) "
        "or TXT file with one SRR/ERR ID per line"
    ),
)
parser.add_argument(
    "-o",
    "--output-dir",
    dest="output_dir",
    default=None,
    help=(
        "Directory to write sample_list.txt. If omitted, writes to "
        "POSEIDON/<Sheet>/<CancerType>/sample_list.txt"
    ),
)
args = parser.parse_args()

# Check if running from correct directory
check_working_directory()

# Determine input type and process accordingly
# Use absolute path to be robust when wrappers change directories
input_path = os.path.abspath(args.metadata_file)
input_ext = os.path.splitext(input_path)[1].lower()
cancer_type = infer_cancer_type_from_filename(input_path)
output_dir_override = os.path.abspath(args.output_dir) if args.output_dir else None

def _write_sample_list(
    mapping: dict,
    sheet_name: str,
    cancer_type: str,
    output_dir_override: Optional[str] = None,
) -> str:
    """Write mapping to sample_list.txt and return the output path.

    mapping format: { sample_id: {"fastq_read1": [...], "fastq_read2": [...]} }
    """
    output_dir = output_dir_override or os.path.join(POSEIDON_DIR, sheet_name, cancer_type)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "sample_list.txt")
    with open(output_file, "w") as fh:
        for sample_id, info in mapping.items():
            fastq_read1_joined = ",".join(info["fastq_read1"]) if info["fastq_read1"] else ""
            # If no R2 entries, write 'NA' for single-end consistency
            fastq_read2_joined = ",".join(info["fastq_read2"]) if info["fastq_read2"] else "NA"
            fh.write(f"{sample_id}\t{fastq_read1_joined}\t{fastq_read2_joined}\n")
    return output_file

def _resolve_biosample_for_srr(srr_id: str) -> Optional[str]:
    """Return BioSample accession (e.g., SAMNxxxxx) for an SRR using sra-tools; fallback to ENA.

    Strategy:
    1) Try `vdb-dump` (sra-tools) in info/JSON modes and extract SAMN[0-9]+.
    2) Fallback to ENA filereport API for biosample_accession.
    """
    # Prefer ENA first for speed; then NCBI efetch XML; then vdb-dump
    try:
        url = (
            "https://www.ebi.ac.uk/ena/portal/api/filereport?"
            f"accession={srr_id}&result=read_run&fields=run_accession,biosample_accession&format=tsv"
        )
        with urlrequest.urlopen(url, timeout=5) as resp:
            body = resp.read().decode("utf-8", errors="ignore")
        lines = [ln.strip() for ln in body.splitlines() if ln.strip()]
        if len(lines) >= 2:
            header = lines[0].split("\t")
            row = lines[1].split("\t")
            if "biosample_accession" in header:
                idx = header.index("biosample_accession")
                biosample = row[idx].strip()
                if biosample:
                    return biosample
    except (urlerror.URLError, TimeoutError, Exception):
        pass

    # NCBI efetch (SRA XML) fallback
    try:
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
            f"db=sra&id={srr_id}&rettype=full&retmode=xml"
        )
        with urlrequest.urlopen(url, timeout=6) as resp:
            xml_text = resp.read().decode("utf-8", errors="ignore")
        m = re.search(r"SAMN\d+", xml_text)
        if m:
            return m.group(0)
    except (urlerror.URLError, TimeoutError, Exception):
        pass

    # vdb-dump in JSON mode
    try:
        proc = subprocess.run(
            ["vdb-dump", srr_id, "-J"], capture_output=True, text=True, check=False, timeout=8
        )
        text = (proc.stdout or "") + (proc.stderr or "")
        m = re.search(r"SAMN\d+", text)
        if m:
            return m.group(0)
    except subprocess.TimeoutExpired:
        pass
    except Exception:
        pass

    # vdb-dump --info mode
    try:
        proc = subprocess.run(
            ["vdb-dump", srr_id, "--info"], capture_output=True, text=True, check=False, timeout=8
        )
        text = (proc.stdout or "") + (proc.stderr or "")
        m = re.search(r"SAMN\d+", text)
        if m:
            return m.group(0)
    except subprocess.TimeoutExpired:
        pass
    except Exception:
        pass

    return None


def process_srr_txt(file_path: str, cancer_type: str, output_dir_override: Optional[str] = None) -> None:
    """Process a TXT file containing one SRR/ERR accession per line.

    Generates POSEIDON/<sheet_name>/<cancer_type>/sample_list.txt where sheet_name is 'SRR'.
    Each SRR becomes its own sample_id row with R1/R2 FASTQ filenames.
    """
    srr_ids = []
    with open(file_path, "r") as fh:
        for raw in fh:
            val = (raw or "").strip()
            if not val:
                continue
            token = val.split()[0]  # allow whitespace-separated values, take first token
            if token.upper().startswith(("SRR", "ERR")) and token[3:].isdigit():
                srr_ids.append(token)
            else:
                # silently skip lines that do not look like SRR/ERR accessions
                continue

    if not srr_ids:
        print(f"No SRR/ERR accessions found in: {file_path}")
        return

    print(f"Resolving BioSample accessions for {len(srr_ids)} SRRs...", flush=True)
    print("Press Ctrl-C to stop early; partial output will be written.", flush=True)
    # Group SRRs by BioSample when available; fallback to SRR as its own sample_id
    mapping: Dict[str, Dict[str, List[str]]] = {}
    layout_cache: Dict[str, str] = {}
    try:
        for idx, srr in enumerate(srr_ids, start=1):
            print(f"  [{idx}/{len(srr_ids)}] {srr}", flush=True)
            biosample = _resolve_biosample_for_srr(srr)
            sample_id = biosample if biosample else srr
            entry = mapping.setdefault(sample_id, {"fastq_read1": [], "fastq_read2": []})
            layout = _resolve_layout_for_run(srr, layout_cache)
            r1 = f"{srr}_1.fastq.gz"
            if r1 not in entry["fastq_read1"]:
                entry["fastq_read1"].append(r1)
            if layout == "PAIRED":
                r2 = f"{srr}_2.fastq.gz"
                if r2 not in entry["fastq_read2"]:
                    entry["fastq_read2"].append(r2)
    except KeyboardInterrupt:
        print("Interrupted by user. Writing partial sample_list.txt...", flush=True)

    out_path = _write_sample_list(
        mapping,
        sheet_name="SRR",
        cancer_type=cancer_type,
        output_dir_override=output_dir_override,
    )

    for sample_id, info in mapping.items():
        # Derive SRR IDs back from fastq names for display
        srrs = sorted({p.split("_")[0] for p in info["fastq_read1"]})
        print(f"Sample: {sample_id} -> Associated SRR IDs: {', '.join(srrs)}")
    print(f"Wrote {len(mapping)} rows to: {out_path}")

# --- Helpers to resolve library layout ---
def _resolve_layout_for_run(run_id: str, cache: Dict[str, str]) -> str:
    """Return 'PAIRED' or 'SINGLE' for a run accession using ENA/Entrez. Caches results."""
    if run_id in cache:
        return cache[run_id]
    layout = _query_layout_via_ena(run_id)
    if not layout:
        layout = _query_layout_via_entrez(run_id)
    if layout not in {"PAIRED", "SINGLE"}:
        layout = "PAIRED"
    cache[run_id] = layout
    return layout

def _query_layout_via_ena(run_id: str) -> Optional[str]:
    try:
        url = (
            "https://www.ebi.ac.uk/ena/portal/api/filereport?"
            f"accession={run_id}&result=read_run&fields=run_accession,library_layout&format=tsv"
        )
        with urlrequest.urlopen(url, timeout=6) as resp:
            body = resp.read().decode("utf-8", errors="ignore")
        lines = [ln.strip() for ln in body.splitlines() if ln.strip()]
        if len(lines) >= 2:
            header = lines[0].split("\t")
            row = lines[1].split("\t")
            if "library_layout" in header:
                idx = header.index("library_layout")
                val = row[idx].strip().upper()
                if val:
                    return val
    except (urlerror.URLError, TimeoutError, Exception):
        return None
    return None

def _query_layout_via_entrez(run_id: str) -> Optional[str]:
    try:
        # esearch+efetch runinfo returns CSV with LibraryLayout column
        cmd = f'esearch -db sra -query "{run_id}" | efetch -format runinfo'
        proc = subprocess.run(cmd, shell=True, capture_output=True, text=True, timeout=10)
        if proc.returncode != 0 or not proc.stdout:
            return None
        import csv as _csv
        for row in _csv.DictReader(proc.stdout.splitlines()):
            if (row.get("Run") or "").strip() == run_id:
                val = (row.get("LibraryLayout") or "").strip().upper()
                if val:
                    return val
        # Fallback: take first row
        rows = list(_csv.DictReader(proc.stdout.splitlines()))
        if rows:
            val = (rows[0].get("LibraryLayout") or "").strip().upper()
            if val:
                return val
    except Exception:
        return None
    return None

if input_ext in {".xlsx", ".xls"}:
    # Expected sheet names to iterate
    target_sheets = ["Tumors", "Controls", "Bulk_CellTypes", "Premalignant"]
    try:
        excel_file = pd.ExcelFile(input_path)
    except Exception as e:
        print(f"Failed to open Excel file: {e}")
        sys.exit(1)

    available = set(excel_file.sheet_names)
    
    # Report sheet status
    print(f"Excel file: {input_path}")
    print(f"Available sheets: {sorted(available)}")
    print(f"Target sheets: {target_sheets}")
    
    found_sheets = []
    missing_sheets = []
    extra_sheets = []
    
    for sheet in target_sheets:
        if sheet in available:
            found_sheets.append(sheet)
        else:
            missing_sheets.append(sheet)
    
    for sheet in available:
        if sheet not in target_sheets:
            extra_sheets.append(sheet)
    
    print(f"\nSheet Analysis:")
    print(f"? Found target sheets: {found_sheets}")
    if missing_sheets:
        print(f"? Missing target sheets: {missing_sheets}")
    if extra_sheets:
        print(f"? Extra sheets (not processed): {extra_sheets}")
    print()
    
    # Process found sheets
    for sheet in found_sheets:
        print(f"Processing sheet: {sheet}")
        df_sheet = pd.read_excel(excel_file, sheet_name=sheet)
        process_dataframe(df_sheet, sheet, cancer_type)
        print()
elif input_ext in {".csv"}:
    # Backwards compatibility: CSV processing writes under POSEIDON/CSV/<CancerType>
    df = pd.read_csv(input_path, delimiter=",")
    # Use a neutral placeholder for sheet to avoid changing CSV behavior
    process_dataframe(df, sheet_name="CSV", cancer_type=cancer_type)
elif input_ext in {".txt"}:
    # New: SRR/ERR list mode (one accession per line)
    process_srr_txt(input_path, cancer_type=cancer_type, output_dir_override=output_dir_override)
else:
    print(f"Unsupported input extension '{input_ext}'. Expected one of: .xlsx, .xls, .csv, .txt")
    sys.exit(1)
