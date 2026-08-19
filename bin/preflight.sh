#!/usr/bin/env bash
# Validate Hellbender, environment, project paths, and Snakemake parsing.

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/.." && pwd)
config_file=
create_envs=false
workflow_prefix=${NGS_WORKFLOW_ENV:-"$HOME/data/conda/envs/hellbender-ngs-workflow"}

usage() {
  printf 'Usage: %s [--create-envs] CONFIG_YAML\n' "$0"
}

while (($#)); do
  case "$1" in
    --create-envs)
      create_envs=true
      shift
      ;;
    -h|--help)
      usage
      exit 0
      ;;
    -*)
      printf 'Unknown option: %s\n' "$1" >&2
      usage >&2
      exit 2
      ;;
    *)
      if [[ -n $config_file ]]; then
        printf 'ERROR: Only one config file may be supplied.\n' >&2
        exit 2
      fi
      config_file=$1
      shift
      ;;
  esac
done

if [[ -z $config_file || ! -r $config_file ]]; then
  printf 'ERROR: Config is not readable: %s\n' "$config_file" >&2
  exit 2
fi
config_file=$(cd -- "$(dirname -- "$config_file")" && pwd)/$(basename -- "$config_file")

if [[ -z ${SLURM_JOB_ID:-} ]] && [[ $(hostname -s) == hellbender-* ]]; then
  printf 'ERROR: Run preflight inside an interactive allocation.\n' >&2
  printf 'Example: srun -p interactive --time=01:00:00 --mem=8G --cpus-per-task=2 --pty bash\n' >&2
  exit 2
fi

if command -v module >/dev/null 2>&1; then
  module purge
  module load miniconda3
fi
if ! command -v conda >/dev/null 2>&1; then
  printf 'ERROR: conda is unavailable.\n' >&2
  exit 2
fi
if [[ ! -d $workflow_prefix/conda-meta ]]; then
  printf 'ERROR: Workflow environment is missing: %s\n' "$workflow_prefix" >&2
  printf 'Run bin/bootstrap_hellbender.sh first.\n' >&2
  exit 2
fi

conda_base=$(conda info --base)
# The module determines this path at runtime.
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
conda activate "$workflow_prefix"

printf 'Validating project configuration...\n'
python "$repo_dir/workflow/lib/pipeline_config.py" --config "$config_file"

if command -v scontrol >/dev/null 2>&1; then
  if ! scontrol show partition general >/dev/null 2>&1; then
    printf 'ERROR: Hellbender partition general is not visible. Check current site policy.\n' >&2
    exit 2
  fi
  printf 'Hellbender general partition visible.\n'
fi

if [[ -d /local/scratch ]]; then
  scratch_parent=/local/scratch/$USER
  mkdir -p "$scratch_parent"
  scratch_test=$(mktemp -d "$scratch_parent/ngs-preflight.XXXXXX")
  rmdir "$scratch_test"
  printf 'Node-local scratch writable: %s\n' "$scratch_parent"
elif [[ $(hostname -s) == hellbender-* ]]; then
  printf 'ERROR: /local/scratch is not available on this compute node.\n' >&2
  exit 2
fi

printf 'Parsing Snakemake workflow...\n'
cd "$repo_dir"
snakemake \
  --snakefile "$repo_dir/workflow/Snakefile" \
  --configfile "$config_file" \
  --profile "$repo_dir/profiles/hellbender-single-node" \
  --list-rules >/dev/null

if [[ $create_envs == true ]]; then
  printf 'Creating rule environments in the shared Snakemake cache...\n'
  snakemake \
    --snakefile "$repo_dir/workflow/Snakefile" \
    --configfile "$config_file" \
    --profile "$repo_dir/profiles/hellbender-single-node" \
    --conda-create-envs-only
fi

printf 'Preflight passed. No analysis jobs were run.\n'
