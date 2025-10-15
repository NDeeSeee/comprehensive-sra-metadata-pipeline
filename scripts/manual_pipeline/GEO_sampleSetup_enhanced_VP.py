import os
import sys
import argparse
from collections import defaultdict

import pandas as pd


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


def _resolve_column_name(df: pd.DataFrame, candidates: list[str]) -> str | None:
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

    for _, row in df.iterrows():
        biosample_raw = row.get(biosample_column)
        run_raw = row.get(run_column)

        if pd.isna(biosample_raw) or pd.isna(run_raw):
            continue

        biosample_id = str(biosample_raw).strip()
        srr_id = str(run_raw).strip()
        if not biosample_id or not srr_id:
            continue

        fastq_r1 = f"{srr_id}_1.fastq.gz"
        fastq_r2 = f"{srr_id}_2.fastq.gz"
        sra_file = f"{srr_id}.sra"

        # Include all rows regardless of file existence
        mapping[biosample_id]["srr_ids"].append(srr_id)
        mapping[biosample_id]["fastq_read1"].append(fastq_r1)
        mapping[biosample_id]["fastq_read2"].append(fastq_r2)

    # Prepare output directory and write results
    output_dir = os.path.join("POSEIDON", sheet_name, cancer_type)
    os.makedirs(output_dir, exist_ok=True)
    output_file = os.path.join(output_dir, "sample_list.txt")

    with open(output_file, "w") as fh:
        for biosample_id, info in mapping.items():
            fastq_read1_joined = ",".join(info["fastq_read1"])
            fastq_read2_joined = ",".join(info["fastq_read2"])
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
    description="Process metadata (CSV or Excel) and organize files by BioSample."
)
parser.add_argument(
    "metadata_file", help="Path to the metadata file (.csv, .xlsx, or .xls)"
)
args = parser.parse_args()

# Determine input type and process accordingly
input_ext = os.path.splitext(args.metadata_file)[1].lower()
cancer_type = infer_cancer_type_from_filename(args.metadata_file)

if input_ext in {".xlsx", ".xls"}:
    # Expected sheet names to iterate
    target_sheets = ["Tumors", "Controls", "Bulk_CellTypes", "Premalignant"]
    try:
        excel_file = pd.ExcelFile(args.metadata_file)
    except Exception as e:
        print(f"Failed to open Excel file: {e}")
        sys.exit(1)

    available = set(excel_file.sheet_names)
    
    # Report sheet status
    print(f"Excel file: {args.metadata_file}")
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
    print(f"✓ Found target sheets: {found_sheets}")
    if missing_sheets:
        print(f"✗ Missing target sheets: {missing_sheets}")
    if extra_sheets:
        print(f"ℹ Extra sheets (not processed): {extra_sheets}")
    print()
    
    # Process found sheets
    for sheet in found_sheets:
        print(f"Processing sheet: {sheet}")
        df_sheet = pd.read_excel(excel_file, sheet_name=sheet)
        process_dataframe(df_sheet, sheet, cancer_type)
        print()
else:
    # Backwards compatibility: CSV processing writes under POSEIDON/CSV/<CancerType>
    df = pd.read_csv(args.metadata_file, delimiter=",")
    # Use a neutral placeholder for sheet to avoid changing CSV behavior
    process_dataframe(df, sheet_name="CSV", cancer_type=cancer_type)
