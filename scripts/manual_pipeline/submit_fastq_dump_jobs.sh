#!/bin/bash
# submit_fastq_dump_jobs.sh  — LSF job shim for SRA→FASTQ
# Uses fasterq-dump + pigz, never deletes .sra/.sralite (Python owns cleanup)

set -euo pipefail

if [ $# -eq 0 ]; then
  echo "Usage: $0 <sra_file_or_sample_list>"
  exit 1
fi

INPUT="$1"
DIR="$(pwd)"

# Helper: submit one SRR
submit_one() {
  local SRR_ID="$1"
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
#    --split-files ensures _1/_2; --threads matches -n; --temp isolates temp I/O
fasterq-dump "\${SRR_ID}" \
  --split-files \
  --threads 8 \
  --temp "\$SCRATCH" \
  --outdir "\$SCRATCH" \
  --progress

# 2) Gzip in parallel (pigz). If pigz missing, fall back to gzip.
if command -v pigz >/dev/null; then
  pigz -p 8 "\$SCRATCH/\${SRR_ID}_1.fastq" 2>/dev/null || true
  [ -f "\$SCRATCH/\${SRR_ID}_2.fastq" ] && pigz -p 8 "\$SCRATCH/\${SRR_ID}_2.fastq" 2>/dev/null || true
else
  gzip "\$SCRATCH/\${SRR_ID}_1.fastq"
  [ -f "\$SCRATCH/\${SRR_ID}_2.fastq" ] && gzip "\$SCRATCH/\${SRR_ID}_2.fastq" || true
fi

# 3) Atomically move finished files into DIR
#    (move then fsync pattern to avoid half-written visibility)
for f in "\$SCRATCH/\${SRR_ID}_1.fastq.gz" "\$SCRATCH/\${SRR_ID}_2.fastq.gz"; do
  if [ -f "\$f" ]; then
    base="\$(basename "\$f")"
    tmp="\$DIR/\${base}.part"
    final="\$DIR/\${base}"
    cp -f "\$f" "\$tmp"
    sync
    mv -f "\$tmp" "\$final"
  fi
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
    # Skip if R1 exists and is non-empty; if paired and R2 exists non-empty, also skip
    if [ -s "${SRR}_1.fastq.gz" ] && { [ ! -f "${SRR}_2.fastq.gz" ] || [ -s "${SRR}_2.fastq.gz" ]; }; then
      echo "✓ ${SRR} FASTQs present; skipping"
      continue
    fi
    submit_one "$SRR" && count=$((count+1))
  done
  echo "Submitted $count jobs."

else
  # Single file or accession
  SAMPLE="$INPUT"
  # If a file path, strip extension to accession
  if [ -f "$SAMPLE" ]; then
    SAMPLE="$(basename "$SAMPLE")"
    SAMPLE="${SAMPLE%.sra}" 
    SAMPLE="${SAMPLE%.sralite}"
  fi
  submit_one "$SAMPLE"
fi