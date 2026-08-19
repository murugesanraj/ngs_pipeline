#!/usr/bin/env bash
# Submit optional, licensed Illumina BCL Convert preprocessing.

set -euo pipefail
umask 0027

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/.." && pwd)
run_dir=
sample_sheet=
output_dir=
partition=general
account=
module_name=
command_name=bcl-convert

usage() {
  cat <<EOF
Usage: $0 --run-dir DIR --sample-sheet FILE --output-dir DIR [options]

Options:
  --partition NAME       Driver partition (default: general)
  --account NAME         Optional investor account
  --module NAME          Optional site module providing BCL Convert
  --command PATH         BCL Convert executable (default: bcl-convert)
EOF
}

while (($#)); do
  case "$1" in
    --run-dir) run_dir=${2:?}; shift 2 ;;
    --sample-sheet) sample_sheet=${2:?}; shift 2 ;;
    --output-dir) output_dir=${2:?}; shift 2 ;;
    --partition) partition=${2:?}; shift 2 ;;
    --account) account=${2:?}; shift 2 ;;
    --module) module_name=${2:?}; shift 2 ;;
    --command) command_name=${2:?}; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ ! -d $run_dir || ! -r $sample_sheet || -z $output_dir ]]; then
  printf 'ERROR: A readable run directory, sample sheet, and output directory are required.\n' >&2
  exit 2
fi
if ! command -v sbatch >/dev/null 2>&1; then
  printf 'ERROR: sbatch is unavailable.\n' >&2
  exit 2
fi

mkdir -p "$output_dir/logs"
run_dir=$(cd -- "$run_dir" && pwd)
sample_sheet=$(cd -- "$(dirname -- "$sample_sheet")" && pwd)/$(basename -- "$sample_sheet")
output_dir=$(cd -- "$output_dir" && pwd)

for value in "$run_dir" "$sample_sheet" "$output_dir" "$module_name" "$command_name"; do
  if [[ $value == *,* || $value == *$'\n'* ]]; then
    printf 'ERROR: Values cannot contain commas or newlines.\n' >&2
    exit 2
  fi
done

args=(
  --parsable
  --job-name=bcl-convert
  --partition="$partition"
  --time=12:00:00
  --cpus-per-task=32
  --mem=64G
  --output="$output_dir/logs/bcl-convert-%j.out"
  --export="ALL,BCL_RUN_DIR=$run_dir,BCL_SAMPLE_SHEET=$sample_sheet,BCL_OUTPUT_DIR=$output_dir,BCL_MODULE=$module_name,BCL_COMMAND=$command_name"
)
if [[ -n $account ]]; then
  args+=(--account="$account")
fi

job_id=$(sbatch "${args[@]}" "$repo_dir/slurm/bcl_convert.sbatch")
printf 'Submitted BCL Convert job %s\n' "$job_id"

