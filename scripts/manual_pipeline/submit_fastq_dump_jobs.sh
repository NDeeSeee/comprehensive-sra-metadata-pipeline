#!/bin/bash
# submit_fastq_dump_jobs.sh  — LSF job shim for SRA→FASTQ
# Uses fasterq-dump + pigz, never deletes .sra/.sralite (Python owns cleanup)

set -euo pipefail

# Normalize and sanitize an accession or filename stem:
# - Convert NBSP (U+00A0) to regular space
# - Trim leading/trailing whitespace
# - Keep only letters and digits (valid for SRR/ERR IDs)
sanitize_id() {
  local s
  s="${1:-}"
  s=$(printf '%s' "$s" | tr '\302\240' ' ')
  s=$(printf '%s' "$s" | sed -E 's/^[[:space:]]+|[[:space:]]+$//g')
  s=$(printf '%s' "$s" | tr -cd '[:alnum:]')
  printf '%s' "$s"
}

if [ $# -eq 0 ]; then
  echo "Usage: $0 <sra_file_or_sample_list>"
  exit 1
fi

INPUT="$1"
DIR="$(pwd)"

# Helper: submit one SRR
submit_one() {
  local SRR_ID_RAW="$1"
  local SRR_ID
  SRR_ID="$(sanitize_id "$SRR_ID_RAW")"
  mkdir -p "$DIR/logs"

  bsub <<EOF
#BSUB -L /bin/bash
#BSUB -W 72:00
#BSUB -n 8
#BSUB -R "span[hosts=1]"
#BSUB -R "rusage[mem=32000]"
#BSUB -M 32000
#BSUB -J fastq_${SRR_ID}
#BSUB -o "$DIR/logs/fastq_${SRR_ID}.out.txt"
#BSUB -e "$DIR/logs/fastq_${SRR_ID}.err.txt"

set -euo pipefail

SRR_ID="${SRR_ID}"
DIR="$DIR"

module load sratoolkit/3.0.10 2>/dev/null || module load sratoolkit/2.11.0 2>/dev/null || true
command -v fasterq-dump >/dev/null || { echo "fasterq-dump not found" >&2; exit 127; }

# Local scratch for temp files (prefer node-local; fall back to /tmp)
SCRATCH="\${LSB_JOB_TMPDIR:-/tmp}/sra_\${SRR_ID}"
mkdir -p "\$SCRATCH"
ulimit -n 4096 || true

cd "\$DIR"

# 1) Convert to uncompressed FASTQ into scratch
#    Prefer local .sra/.sralite if present; else use accession
SRC="\${SRR_ID}"
if [ -f "\$DIR/\${SRR_ID}.sra" ]; then
  SRC="\$DIR/\${SRR_ID}.sra"
elif [ -f "\$DIR/\${SRR_ID}.sralite" ]; then
  SRC="\$DIR/\${SRR_ID}.sralite"
fi

#    --split-files ensures _1/_2; --threads matches -n; --temp isolates temp I/O
fasterq-dump "\${SRC}" \
  --split-files \
  --threads 8 \
  --temp "\$SCRATCH" \
  --outdir "\$SCRATCH" \
  --progress

# 2) Gzip in parallel (pigz). Handle any files that start with SRR_ID
shopt -s nullglob
FASTQS=( "\$SCRATCH/\${SRR_ID}"*.fastq )
if command -v pigz >/dev/null; then
  if [ \${#FASTQS[@]} -gt 0 ]; then
    pigz -p 8 "\${FASTQS[@]}" 2>/dev/null || true
  fi
else
  for fq in "\${FASTQS[@]}"; do
    gzip "\$fq" || true
  done
fi

# 3) Atomically move finished files into DIR
#    (move then fsync pattern to avoid half-written visibility)
for f in "\$SCRATCH/\${SRR_ID}"*.fastq.gz; do
  [ -f "\$f" ] || continue
  base="\$(basename "\$f")"
  # Standardize single-end name if toolkit produced \${SRR_ID}.fastq.gz
  if [ "\$base" = "\${SRR_ID}.fastq.gz" ]; then
    base="\${SRR_ID}_1.fastq.gz"
  fi
  tmp="\$DIR/\${base}.part"
  final="\$DIR/\${base}"
  cp -f "\$f" "\$tmp"
  sync
  mv -f "\$tmp" "\$final"
done

# 4) Leave .sra/.sralite in place — Python validates & deletes safely.
#    Success exit code signals Python to perform thorough validation.
EOF
}

# Branch: sample_list.txt
if [ "$(basename "$INPUT")" = "sample_list.txt" ]; then
  echo "Processing sample_list.txt..."
  if [ ! -f "$INPUT" ]; then
    echo "ERROR: sample_list.txt not found" >&2
    exit 1
  fi

  # pull SRR/ERR from both R1 and R2 columns, tolerate commas
  SRR_IDS=$(
    awk -F'\t' 'NF>=3 {print $2","$3}' "$INPUT" \
    | tr ',' '\n' \
    | sed -E 's/(_1|_2)?\.fastq(\.gz)?$//' \
    | grep -E '^(SRR|ERR)[0-9]+' \
    | sort -u
  )

  count=0
  for SRR in $SRR_IDS; do
    SRR_CLEAN="$(sanitize_id "$SRR")"
    # Skip if R1 exists and is non-empty; if paired and R2 exists non-empty, also skip
    if [ -s "${SRR_CLEAN}_1.fastq.gz" ] && { [ ! -f "${SRR_CLEAN}_2.fastq.gz" ] || [ -s "${SRR_CLEAN}_2.fastq.gz" ]; }; then
      echo "✓ ${SRR_CLEAN} FASTQs present; skipping"
      continue
    fi
    submit_one "$SRR_CLEAN" && count=$((count+1))
  done
  echo "Submitted $count jobs."

else
  # Single file or accession
  SAMPLE_RAW="$INPUT"
  # If a file path, strip extension to accession
  if [ -f "$SAMPLE_RAW" ]; then
    SAMPLE_RAW="$(basename "$SAMPLE_RAW")"
    SAMPLE_RAW="${SAMPLE_RAW%.sra}"
    SAMPLE_RAW="${SAMPLE_RAW%.sralite}"
  fi
  SAMPLE="$(sanitize_id "$SAMPLE_RAW")"
  submit_one "$SAMPLE"
fi