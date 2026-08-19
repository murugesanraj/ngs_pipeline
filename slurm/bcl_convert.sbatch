#!/usr/bin/env bash
#SBATCH --job-name=bcl-convert
#SBATCH --partition=general
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --export=ALL

set -euo pipefail
umask 0027

: "${BCL_RUN_DIR:?BCL_RUN_DIR is required}"
: "${BCL_SAMPLE_SHEET:?BCL_SAMPLE_SHEET is required}"
: "${BCL_OUTPUT_DIR:?BCL_OUTPUT_DIR is required}"
BCL_COMMAND=${BCL_COMMAND:-bcl-convert}

module purge
if [[ -n ${BCL_MODULE:-} ]]; then
  module load "$BCL_MODULE"
fi
if ! command -v "$BCL_COMMAND" >/dev/null 2>&1; then
  printf 'ERROR: BCL Convert executable is unavailable: %s\n' "$BCL_COMMAND" >&2
  printf 'Use --module or --command with an institutionally approved installation.\n' >&2
  exit 127
fi

mkdir -p "$BCL_OUTPUT_DIR"
"$BCL_COMMAND" \
  --bcl-input-directory "$BCL_RUN_DIR" \
  --sample-sheet "$BCL_SAMPLE_SHEET" \
  --output-directory "$BCL_OUTPUT_DIR" \
  --force

find "$BCL_OUTPUT_DIR" -maxdepth 1 -type f -name '*.fastq.gz' -print0 \
  | sort -z \
  | xargs -0 --no-run-if-empty sha256sum \
  > "$BCL_OUTPUT_DIR/fastq.sha256"

if [[ ! -s "$BCL_OUTPUT_DIR/fastq.sha256" ]]; then
  printf 'ERROR: BCL Convert produced no top-level FASTQ files in %s\n' "$BCL_OUTPUT_DIR" >&2
  exit 1
fi
