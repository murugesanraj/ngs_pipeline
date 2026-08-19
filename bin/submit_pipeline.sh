#!/usr/bin/env bash
# Submit a lightweight Snakemake controller job; analysis rules submit via SLURM.

set -euo pipefail
umask 0027

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/.." && pwd)
config_file=
profile_dir="$repo_dir/profiles/hellbender"
driver_partition=general
driver_account=
workflow_prefix=${NGS_WORKFLOW_ENV:-"$HOME/data/conda/envs/hellbender-ngs-workflow"}
single_node=false

usage() {
  cat <<EOF
Usage: $0 --config FILE [options]

Options:
  --profile DIR              Snakemake profile (default: profiles/hellbender)
  --driver-partition NAME    Partition for controller allocation (default: general)
  --driver-account NAME      Optional investor account for controller allocation
  --workflow-env PATH        Workflow Conda prefix
  --single-node              Run all rules locally inside one 32-CPU/240-GB allocation
  -h, --help                 Show this help
EOF
}

while (($#)); do
  case "$1" in
    --config)
      config_file=${2:?Missing value for --config}
      shift 2
      ;;
    --profile)
      profile_dir=${2:?Missing value for --profile}
      shift 2
      ;;
    --driver-partition)
      driver_partition=${2:?Missing value for --driver-partition}
      shift 2
      ;;
    --driver-account)
      driver_account=${2:?Missing value for --driver-account}
      shift 2
      ;;
    --workflow-env)
      workflow_prefix=${2:?Missing value for --workflow-env}
      shift 2
      ;;
    --single-node)
      single_node=true
      profile_dir="$repo_dir/profiles/hellbender-single-node"
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    *)
      printf 'Unknown option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
  esac
done

if [[ -z $config_file || ! -r $config_file ]]; then
  printf 'ERROR: --config must name a readable YAML file.\n' >&2
  exit 2
fi
if [[ ! -d $profile_dir ]]; then
  printf 'ERROR: Profile directory does not exist: %s\n' "$profile_dir" >&2
  exit 2
fi
if [[ ! -d $workflow_prefix/conda-meta ]]; then
  printf 'ERROR: Workflow environment is missing: %s\n' "$workflow_prefix" >&2
  exit 2
fi
if ! command -v sbatch >/dev/null 2>&1; then
  printf 'ERROR: sbatch is unavailable; submit this command on Hellbender.\n' >&2
  exit 2
fi

for value in "$config_file" "$profile_dir" "$workflow_prefix" "$repo_dir"; do
  if [[ $value == *,* || $value == *$'\n'* ]]; then
    printf 'ERROR: Exported paths cannot contain commas or newlines: %s\n' "$value" >&2
    exit 2
  fi
done

config_file=$(cd -- "$(dirname -- "$config_file")" && pwd)/$(basename -- "$config_file")
profile_dir=$(cd -- "$profile_dir" && pwd)
workflow_prefix=$(cd -- "$workflow_prefix" && pwd)
mkdir -p "$repo_dir/logs/driver"

sbatch_args=(
  --parsable
  --job-name=ngs-driver
  --partition="$driver_partition"
  --time=2-00:00:00
  --output="$repo_dir/logs/driver/%x-%j.out"
  --export="ALL,NGS_CONFIG=$config_file,NGS_PROFILE=$profile_dir,NGS_WORKFLOW_ENV=$workflow_prefix,NGS_REPO_DIR=$repo_dir"
)

if [[ -n $driver_account ]]; then
  sbatch_args+=(--account="$driver_account")
fi
if [[ $single_node == true ]]; then
  sbatch_args+=(--cpus-per-task=32 --mem=240G)
else
  sbatch_args+=(--cpus-per-task=2 --mem=8G)
fi

job_id=$(sbatch "${sbatch_args[@]}" "$repo_dir/slurm/driver.sbatch")
printf 'Submitted workflow controller job %s\n' "$job_id"
printf 'Monitor with: squeue -j %s\n' "$job_id"

