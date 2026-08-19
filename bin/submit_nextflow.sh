#!/usr/bin/env bash
# Submit a lightweight Nextflow controller allocation on Hellbender.

set -euo pipefail
umask 0027

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/.." && pwd)
config_file=
profile=hellbender
extra_config=
driver_partition=general
driver_account=
nextflow_prefix=${NGS_NEXTFLOW_ENV:-"$HOME/data/conda/envs/hellbender-ngs-nextflow"}
single_node=false

usage() {
  cat <<EOF
Usage: $0 --config FILE [options]

Options:
  --profile NAME             Nextflow profile (default: hellbender)
  --nextflow-config FILE     Optional priority/local Nextflow config
  --driver-partition NAME    Controller partition (default: general)
  --driver-account NAME      Optional investor account for the controller
  --nextflow-env PATH        Nextflow Conda prefix
  --single-node              Run all tasks inside one 32-CPU/240-GB allocation
  -h, --help                 Show this help
EOF
}

while (($#)); do
  case "$1" in
    --config) config_file=${2:?Missing value for --config}; shift 2 ;;
    --profile) profile=${2:?Missing value for --profile}; shift 2 ;;
    --nextflow-config) extra_config=${2:?Missing value for --nextflow-config}; shift 2 ;;
    --driver-partition) driver_partition=${2:?Missing value}; shift 2 ;;
    --driver-account) driver_account=${2:?Missing value}; shift 2 ;;
    --nextflow-env) nextflow_prefix=${2:?Missing value}; shift 2 ;;
    --single-node)
      single_node=true
      profile=hellbender_single_node
      shift
      ;;
    -h|--help) usage; exit 0 ;;
    *) printf 'Unknown option: %s\n' "$1" >&2; usage >&2; exit 2 ;;
  esac
done

if [[ -z $config_file || ! -r $config_file ]]; then
  printf 'ERROR: --config must name a readable YAML file.\n' >&2
  exit 2
fi
if [[ -n $extra_config && ! -r $extra_config ]]; then
  printf 'ERROR: Nextflow config is not readable: %s\n' "$extra_config" >&2
  exit 2
fi
if [[ ! -d $nextflow_prefix/conda-meta ]]; then
  printf 'ERROR: Nextflow environment is missing: %s\n' "$nextflow_prefix" >&2
  exit 2
fi
if ! command -v sbatch >/dev/null 2>&1; then
  printf 'ERROR: sbatch is unavailable; submit this command on Hellbender.\n' >&2
  exit 2
fi

config_file=$(cd -- "$(dirname -- "$config_file")" && pwd)/$(basename -- "$config_file")
nextflow_prefix=$(cd -- "$nextflow_prefix" && pwd)
if [[ -n $extra_config ]]; then
  extra_config=$(cd -- "$(dirname -- "$extra_config")" && pwd)/$(basename -- "$extra_config")
fi
for value in "$config_file" "$extra_config" "$nextflow_prefix" "$repo_dir" "$profile"; do
  if [[ $value == *,* || $value == *$'\n'* ]]; then
    printf 'ERROR: Exported values cannot contain commas or newlines.\n' >&2
    exit 2
  fi
done

mkdir -p "$repo_dir/logs/nextflow-driver"
sbatch_args=(
  --parsable
  --job-name=ngs-nextflow
  --partition="$driver_partition"
  --time=2-00:00:00
  --output="$repo_dir/logs/nextflow-driver/%x-%j.out"
  --export="ALL,NGS_CONFIG=$config_file,NGS_NXF_PROFILE=$profile,NGS_NXF_CONFIG=$extra_config,NGS_NEXTFLOW_ENV=$nextflow_prefix,NGS_REPO_DIR=$repo_dir"
)
if [[ -n $driver_account ]]; then
  sbatch_args+=(--account="$driver_account")
fi
if [[ $single_node == true ]]; then
  sbatch_args+=(--cpus-per-task=32 --mem=240G)
else
  sbatch_args+=(--cpus-per-task=2 --mem=8G)
fi

job_id=$(sbatch "${sbatch_args[@]}" "$repo_dir/slurm/nextflow_driver.sbatch")
printf 'Submitted Nextflow controller job %s\n' "$job_id"
printf 'Monitor with: squeue -j %s\n' "$job_id"
