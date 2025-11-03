# Project-specific aliases and helpers for POSEIDON
# Usage: source this file in your shell session
#   source scripts/project_aliases.sh

# Resolve project root based on this file's location
_project_aliases_resolve_root() {
  local this_dir
  this_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  echo "$(cd "${this_dir}/.." && pwd)"
}

export POSEIDON_PROJECT_ROOT="$(_project_aliases_resolve_root)"
export POSEIDON_MANUAL_DIR="${POSEIDON_PROJECT_ROOT}/scripts/manual_pipeline"
export POSEIDON_DIR="/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON"

# --- Commands ---

sra-summary() {
  # Generate markdown summary of SRA Excel metadata files
  python3 "${POSEIDON_MANUAL_DIR}/generate_sra_summary.py" "$@"
}


srr_to_geo() {
  # Build sample_list.txt under POSEIDON. Conditional copy/merge based on input type.
  # Usage:
  #   srr-to-geo <metadata_or_base_sample_list> [extra_sheet_or_table ...]
  # Behavior:
  #   - Excel (.xlsx/.xls): generate per-sheet sample_list.txt under POSEIDON/<Sheet>/<CancerType>/; do NOT copy to CWD; skip merge.
  #   - CSV/TXT metadata (.csv or SRR list .txt): also copy generated sample_list.txt into CWD; if extras provided, write sample_list.merged.txt in CWD.
  #   - Merge-only mode: if first arg is an existing sample_list.txt, skip generator and merge extras into CWD/sample_list.merged.txt.
  local prev_dir input_arg abs_input
  prev_dir="$(pwd)"
  input_arg="$1"
  shift || true
  # Collect any extra sheets (zero or more)
  local extra_sheets=()
  while [ $# -gt 0 ]; do
    extra_sheets+=("$1")
    shift
  done

  # Resolve first argument to absolute path if provided
  if [ -n "$input_arg" ]; then
    if [ "${input_arg#*/}" != "$input_arg" ]; then
      # contains a slash, resolve relative to current dir
      abs_input="$(cd "$(dirname "$input_arg")" && pwd)/$(basename "$input_arg")"
    else
      # no slash, treat as file in current dir
      abs_input="${prev_dir}/${input_arg}"
    fi
  fi

  # Determine input kind and whether we are in merge-only mode
  local input_kind="unknown" merge_only=0
  if [ -n "$abs_input" ] && [ -f "$abs_input" ]; then
    case "$abs_input" in
      *.xlsx|*.xls) input_kind="excel" ;;
      *.csv)        input_kind="csv" ;;
      *.txt)        input_kind="txt" ;;
      *)            input_kind="unknown" ;;
    esac
    if [ "$(basename "$abs_input")" = "sample_list.txt" ]; then
      merge_only=1
    fi
  fi

  local rc wrote_path
  if [ "${merge_only}" -eq 1 ]; then
    # Skip generator; we'll only perform a merge using the provided base sample_list.txt
    rc=0
    wrote_path=""
  else
    cd "${POSEIDON_DIR}" || { echo "Failed to cd to ${POSEIDON_DIR}" >&2; return 1; }
    # Run generator and capture its stdout to detect where it wrote the file
    local tmp_log
    tmp_log="$(mktemp)"
    if [ -n "$abs_input" ]; then
      python3 "${POSEIDON_PROJECT_ROOT}/scripts/manual_pipeline/generate_geo_sample_lists.py" "$abs_input" | tee "$tmp_log"
      rc=${PIPESTATUS[0]}
    else
      python3 "${POSEIDON_PROJECT_ROOT}/scripts/manual_pipeline/generate_geo_sample_lists.py" | tee "$tmp_log"
      rc=${PIPESTATUS[0]}
    fi
    # Try to parse the output path from the generator message (single-output modes)
    wrote_path="$(awk -F': ' '/^Wrote .* to: /{print $2}' "$tmp_log" | tail -n 1)"
    rm -f "$tmp_log"
  fi

  # Only copy when generator succeeded
  if [ ${rc:-0} -eq 0 ]; then
    # Fallback: pick the most recent sample_list.txt under SRR/ (only relevant for SRR mode)
    if [ -z "$wrote_path" ] || [ ! -f "$wrote_path" ]; then
      wrote_path="$(ls -t "${POSEIDON_DIR}/SRR"/*/sample_list.txt 2>/dev/null | head -n 1)"
    fi
    # Copy into CWD only for CSV/TXT inputs (quick-iteration modes)
    if [ "${merge_only}" -ne 1 ] && { [ "${input_kind}" = "csv" ] || [ "${input_kind}" = "txt" ]; }; then
      if [ -n "$wrote_path" ] && [ -f "$wrote_path" ]; then
        cp -f "$wrote_path" "${prev_dir}/sample_list.txt"
        echo "Copied generated sample list to: ${prev_dir}/sample_list.txt"
      else
        echo "Warning: Could not locate generated sample_list.txt" >&2
      fi
    fi
  fi

  # If extra sheets were provided, merge them into a new file in caller directory
  # Allowed when: merge-only mode OR input was CSV/TXT. Skipped for Excel inputs.
  if [ ${#extra_sheets[@]} -gt 0 ] && { [ "${merge_only}" -eq 1 ] || [ "${input_kind}" = "csv" ] || [ "${input_kind}" = "txt" ]; }; then
    # If generator failed, proceed if caller already has a sample_list.txt
    if [ ${rc:-0} -ne 0 ]; then
      if [ -f "${prev_dir}/sample_list.txt" ]; then
        echo "Warning: generator failed (rc=${rc}); using existing sample_list.txt for merge" >&2
      else
        cd "${prev_dir}" || true
        return ${rc}
      fi
    fi
    # Resolve extras to absolute paths for robust handling
    local resolved_extras=()
    local p
    for p in "${extra_sheets[@]}"; do
      if [ "${p#*/}" != "$p" ]; then
        resolved_extras+=("$(cd "$(dirname "$p")" && pwd)/$(basename "$p")")
      else
        resolved_extras+=("${prev_dir}/$p")
      fi
    done

    # Determine base for merge: merge-only uses provided sample_list.txt; otherwise use the CWD copy
    local base_for_merge out_merge
    if [ "${merge_only}" -eq 1 ]; then
      base_for_merge="$abs_input"
    else
      base_for_merge="${prev_dir}/sample_list.txt"
    fi
    out_merge="${prev_dir}/sample_list.merged.txt"

    # Keep a copy of the raw file alongside
    if [ -f "${base_for_merge}" ]; then
      cp -f "${base_for_merge}" "${prev_dir}/sample_list.raw.txt"
    fi

    # Merge into sample_list.merged.txt (do not overwrite the 3-col file by default)
    python3 - "$base_for_merge" "$out_merge" "${resolved_extras[@]}" <<'PY'
import sys, os, csv, re
from collections import OrderedDict, defaultdict

if len(sys.argv) < 3:
    sys.exit("Usage: merge.py <base_sample_list.tsv> <out.tsv> [extra1 ...]")

base_path, out_path, *extra_paths = sys.argv[1:]

def read_base(path):
    rows = []
    samn_to_idx = {}
    srr_to_samn = {}
    with open(path, 'r', newline='') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 3:
                # pad to 3 columns
                parts = parts + [''] * (3 - len(parts))
            biosample, r1, r2 = parts[0], parts[1], parts[2]
            row = OrderedDict([('BioSample', biosample), ('R1', r1), ('R2', r2)])
            idx = len(rows)
            rows.append(row)
            samn_to_idx[biosample] = idx
            # build SRR/ERR mapping from filenames in R1/R2 columns
            def extract_srrs(v):
                ids = []
                for token in filter(None, [t.strip() for t in v.split(',')]):
                    # Accept SRR/ERR anywhere in the basename (handles *_1.fastq.gz)
                    base = os.path.basename(token)
                    m = re.search(r'(SRR\d+|ERR\d+)', base)
                    if m:
                        ids.append(m.group(1))
                return ids
            for srr in set(extract_srrs(r1) + extract_srrs(r2)):
                srr_to_samn[srr] = biosample
    return rows, samn_to_idx, srr_to_samn

def sniff_delim(sample_line):
    return ',' if sample_line.count(',') > sample_line.count('\t') else '\t'

def read_table(path):
    with open(path, 'r', newline='') as f:
        head = f.readline().rstrip('\n')
        if not head:
            return []
        delim = sniff_delim(head)
        # Detect headerless 3-col sample_list format by checking first line tokens
        tokens = head.split(delim)
        is_sample_list_like = False
        if len(tokens) >= 3:
            t0 = tokens[0].strip()
            t1 = tokens[1].strip()
            t2 = tokens[2].strip()
            if re.match(r'^(SAM[NDE]|SAME|ERS|DRS)\w+', t0) and \
               re.search(r'_(1|R1)\.fastq\.gz$', t1) and \
               re.search(r'_(2|R2)\.fastq\.gz$', t2):
                is_sample_list_like = True
        f.seek(0)
        if is_sample_list_like:
            rows = []
            for line in f:
                s = line.rstrip('\n')
                if not s:
                    continue
                parts = s.split(delim)
                if len(parts) < 3:
                    continue
                rows.append({'BioSample': parts[0].strip(), 'R1': parts[1].strip(), 'R2': parts[2].strip()})
            return rows
        # Fallback: treat as headered CSV/TSV via DictReader
        reader = csv.DictReader(f, delimiter=delim)
        return list(reader)

def normalize_header(h):
    return re.sub(r'[^a-z0-9]+', ' ', h.strip().lower())

def canonicalize_key(key):
    n = normalize_header(key)
    tokens = set(n.split())
    # Canonicalize BioSample column variants
    if 'samn' in tokens or 'biosample' in tokens:
        return 'BioSample'
    # Canonicalize read columns
    if n in {'r1', 'read 1'} or 'r1' in tokens:
        return 'R1'
    if n in {'r2', 'read 2'} or 'r2' in tokens:
        return 'R2'
    return key

def find_cols(headers):
    norm = [normalize_header(h) for h in headers]
    # Candidates for BioSample
    samn_keys = {'biosample', 'biosample accession', 'samn', 'biosample id'}
    biosample_cols = [headers[i] for i, n in enumerate(norm) if any(k in n.split() for k in samn_keys)]
    # Candidates for SRR/ERR or FASTQ filename columns (R1/R2)
    srr_like = []
    for i, n in enumerate(norm):
        tokens = set(n.split())
        if 'srr' in tokens or 'err' in tokens or n in {'r1', 'r2'} or 'fastq' in tokens:
            srr_like.append(headers[i])
    return biosample_cols, srr_like

base_rows, samn_to_idx, srr_to_samn = read_base(base_path)

def ensure_row_for_samn(biosample):
    # For empty biosample, always create a fresh row (do not index by '')
    if biosample and biosample in samn_to_idx:
        return base_rows[samn_to_idx[biosample]]
    idx = len(base_rows)
    row = OrderedDict([('BioSample', biosample), ('R1', ''), ('R2', '')])
    base_rows.append(row)
    if biosample:
        samn_to_idx[biosample] = idx
    return row

def add_field(row, key, value, tag):
    if value is None:
        return
    val = str(value).strip()
    if val == '':
        return
    # Canonicalize incoming keys to core columns when possible
    ckey = canonicalize_key(key)
    if ckey == 'BioSample':
        if not row.get('BioSample'):
            row['BioSample'] = val
        return
    if ckey in {'R1', 'R2'}:
        if not row.get(ckey):
            row[ckey] = val
        return
    if ckey not in row:
        row[ckey] = val
        return
    # If same value or row has empty, fill in
    if not row[ckey]:
        row[ckey] = val
        return
    if row[ckey] == val:
        return
    # Collision: create suffixed column
    suffix = 1
    new_key = f"{ckey}_{tag}"
    while new_key in row:
        suffix += 1
        new_key = f"{ckey}_{tag}{suffix}"
    row[new_key] = val

unmatched_added_total = 0
for idx, extra_path in enumerate(extra_paths, start=1):
    sheet_rows = read_table(extra_path)
    if not sheet_rows:
        continue
    headers = list(sheet_rows[0].keys())
    biosample_cols, srr_cols = find_cols(headers)
    tag = f"sheet{idx}"
    for r in sheet_rows:
        matched_samn = []
        # Try BioSample columns first
        for col in biosample_cols:
            v = (r.get(col) or '').strip()
            if not v:
                continue
            # Extract one or more SAMN tokens from the cell
            for tok in re.split(r'[\s,;|]+', v):
                m = re.search(r'(SAMN\d+)', tok)
                if m:
                    matched_samn.append(m.group(1))
        # If none, try SRR columns
        if not matched_samn:
            srrs = set()
            for col in srr_cols:
                val = r.get(col) or ''
                for token in re.split(r'[\s,;|]+', val.strip()):
                    base = os.path.basename(token)
                    m = re.search(r'(SRR\d+|ERR\d+)', base)
                    if m:
                        srrs.add(m.group(1))
            for s in sorted(srrs):
                bs = srr_to_samn.get(s)
                if bs:
                    matched_samn.append(bs)
        # Fallback: single row with unknown biosample
        if not matched_samn:
            matched_samn = ['']
        # Deduplicate, preserve order
        seen_bs = set()
        uniq_samn = []
        for bs in matched_samn:
            if bs not in seen_bs:
                seen_bs.add(bs)
                uniq_samn.append(bs)
        for bs in uniq_samn:
            if not bs:
                unmatched_added_total += 1
            row = ensure_row_for_samn(bs)
            # Merge all fields from the sheet row
            for k, v in r.items():
                add_field(row, k, v, tag)

# Compose output columns: strictly 3-column output for pipeline compatibility
all_keys = ['BioSample', 'R1', 'R2']

def extract_srrs_from_row(row):
    vals = []
    for key in ('R1', 'R2'):
        v = row.get(key, '') or ''
        if v:
            for token in v.split(','):
                base = os.path.basename(token.strip())
                m = re.search(r'(SRR\d+|ERR\d+)', base)
                if m:
                    vals.append(m.group(1))
    # Stable unique order
    seen = set()
    uniq = []
    for s in vals:
        if s not in seen:
            seen.add(s)
            uniq.append(s)
    return uniq

def _stable_unique_tokens_csv(val: str):
    tokens = []
    seen = set()
    for t in (val or '').split(','):
        t2 = t.strip()
        if not t2:
            continue
        if t2 not in seen:
            seen.add(t2)
            tokens.append(t2)
    return tokens

# Aggregate multi-SRR per BioSample on a single row with comma-separated R1/R2
with open(out_path, 'w', newline='') as out_f:
    writer = csv.writer(out_f, delimiter='\t', lineterminator='\n')
    for row in base_rows:
        biosample = (row.get('BioSample') or '').strip()
        if not biosample:
            continue
        r1_list = _stable_unique_tokens_csv(row.get('R1', '') or '')
        r2_list = _stable_unique_tokens_csv(row.get('R2', '') or '')
        r1 = ','.join(r1_list)
        r2 = ','.join(r2_list)
        writer.writerow([biosample, r1, r2])

print(f"Merged sheet(s) -> {out_path}")
print(f"Unmatched rows added without BioSample: {unmatched_added_total}")
PY
    merge_rc=$?
    if [ ${merge_rc:-0} -eq 0 ]; then
      echo "Wrote merged sample list to: ${prev_dir}/sample_list.merged.txt"
    else
      echo "Error: Merge failed with status ${merge_rc}" >&2
    fi
  fi

  local final_rc=${rc:-0}
  if [ ${merge_rc:-0} -ne 0 ]; then
    final_rc=${merge_rc}
  fi
  cd "${prev_dir}" || true
  return ${final_rc}
}

# Backwards-compatible alias with hyphenated name
alias srr-to-geo=srr_to_geo

fastq-workflow() {
  # Download SRAs, convert to FASTQ, and write status for a cancer directory
  # Usage: fastq-workflow <cancer_directory>
  python3 "${POSEIDON_MANUAL_DIR}/sample_list_to_fastq_workflow.py" "$@"
}

overall-status() {
  # Compute overall pipeline status table
  # Usage: overall-status [root_dir] [out_md]
  python3 "${POSEIDON_MANUAL_DIR}/compute_overall_status.py" "$@"
}

# Convenience wrapper for the parallel batch runner over Controls/Tumors/Premalignant
run_fastq_workflows() {
  python3 "/data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/ValeriiGitRepo/scripts/automated_pipeline/run_fastq_workflows.py" "$@"
}
alias run-fastq-workflows=run_fastq_workflows
alias rwf=run_fastq_workflows

commands_help() {
  cat <<'EOF'
Project commands (source scripts/project_aliases.sh first):

- sra-summary
  Description: Summarize SRA Excel metadata in SRAMetadataFiles into a markdown table.
  Run from: anywhere (wrapper resolves paths)
  Example: sra-summary

- srr-to-geo <metadata_or_base_sample_list> [extra_sheet ...]
  Description:
    - Excel (.xlsx/.xls): write per-sheet sample_list.txt under POSEIDON/<Sheet>/<CancerType>/; no CWD files and no auto-merge.
    - CSV/TXT: also copy generated sample_list.txt to CWD; if extras provided, write sample_list.merged.txt in CWD.
    - Merge-only: if first arg is an existing sample_list.txt, skip generation and merge extras into CWD/sample_list.merged.txt.
  Run from: anywhere; command cd's to /data/salomonis-archive/FASTQs/NCI-R01/POSEIDON as required.
  Examples:
    srr-to-geo Tongue-SRA.xlsx
    srr-to-geo SRRs.txt
    srr-to-geo SRRs.txt sample_sheet.txt
    srr-to-geo /data/.../Controls/<CancerType>/sample_list.txt extra.csv
  HINT: Before running: conda activate sra-metadata

- fastq-workflow <cancer_directory>
  Description: For a cancer directory with sample_list.txt, download SRAs, convert to FASTQ, and write sample_list.with_status.txt
  Run from: anywhere (pass absolute or relative path)
  Example: fastq-workflow /data/salomonis-archive/FASTQs/NCI-R01/POSEIDON/Tumors/Vagina/
  HINT: Before running: conda activate sra-metadata

- overall-status [root_dir] [out_md]
  Description: Compute a markdown table summarizing pipeline status across directories.
  Run from: anywhere (defaults to project ROOT inside script)
  Example: overall-status

- run_fastq_workflows [options]
  Description: Batch-run fastq-workflow across Controls/Tumors/Premalignant in parallel.
  Run from: anywhere
  Example: run_fastq_workflows --dry-run
  Alias: run-fastq-workflows, rwf

Notes:
- To enable these commands in a shell session, run: source scripts/project_aliases.sh
- These are project-scoped helpers; they do not modify your global shell config.
EOF
}
