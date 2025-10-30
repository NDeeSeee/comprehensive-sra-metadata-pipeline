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

srr-to-geo() {
  # Build sample_list.txt files from metadata (CSV/XLSX) under POSEIDON directory
  # Script must run from ${POSEIDON_DIR}; this wrapper ensures that
  local prev_dir
  prev_dir="$(pwd)"
  cd "${POSEIDON_DIR}" || { echo "Failed to cd to ${POSEIDON_DIR}" >&2; return 1; }
  python3 "${POSEIDON_PROJECT_ROOT}/scripts/manual_pipeline/generate_geo_sample_lists.py" "$@"
  local rc=$?
  cd "${prev_dir}" || true
  return ${rc}
}

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

commands_help() {
  cat <<'EOF'
Project commands (source scripts/project_aliases.sh first):

- sra-summary
  Description: Summarize SRA Excel metadata in SRAMetadataFiles into a markdown table.
  Run from: anywhere (wrapper resolves paths)
  Example: sra-summary

- srr-to-geo <metadata_file>
  Description: Create POSEIDON-compatible sample_list.txt from CSV/XLSX metadata.
  Run from: anywhere; command temporarily cd's to /data/salomonis-archive/FASTQs/NCI-R01/POSEIDON as required by the script
  Example: srr-to-geo data/manual_metadata/Pancreas_metadata.xlsx
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

Notes:
- To enable these commands in a shell session, run: source scripts/project_aliases.sh
- These are project-scoped helpers; they do not modify your global shell config.
EOF
}
