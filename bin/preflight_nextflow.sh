#!/usr/bin/env bash
# Validate Hellbender paths, project configuration, and the Nextflow graph.

set -euo pipefail

if (($# != 1)); then
  printf 'Usage: %s CONFIG_YAML\n' "$0" >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/.." && pwd)
config_file=$1
nextflow_prefix=${NGS_NEXTFLOW_ENV:-"$HOME/data/conda/envs/hellbender-ngs-nextflow"}

if [[ ! -r $config_file ]]; then
  printf 'ERROR: Config is not readable: %s\n' "$config_file" >&2
  exit 2
fi
config_file=$(cd -- "$(dirname -- "$config_file")" && pwd)/$(basename -- "$config_file")

if [[ -z ${SLURM_JOB_ID:-} ]] && [[ $(hostname -s) == hellbender-* ]]; then
  printf 'ERROR: Run preflight inside an interactive allocation.\n' >&2
  exit 2
fi
if command -v module >/dev/null 2>&1; then
  module purge
  module load miniconda3
fi
if [[ ! -d $nextflow_prefix/conda-meta ]]; then
  printf 'ERROR: Nextflow environment is missing: %s\n' "$nextflow_prefix" >&2
  printf 'Run bin/bootstrap_nextflow.sh first.\n' >&2
  exit 2
fi

conda_base=$(conda info --base)
# The module determines this path at runtime.
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
conda activate "$nextflow_prefix"

python "$repo_dir/workflow/lib/pipeline_config.py" --config "$config_file"
if command -v scontrol >/dev/null 2>&1; then
  scontrol show partition general >/dev/null
  printf 'Hellbender general partition visible.\n'
fi
if [[ -d /local/scratch ]]; then
  scratch_parent=/local/scratch/$USER
  mkdir -p "$scratch_parent"
  scratch_test=$(mktemp -d "$scratch_parent/ngs-nextflow-preflight.XXXXXX")
  rmdir "$scratch_test"
  printf 'Node-local scratch writable: %s\n' "$scratch_parent"
fi

export NXF_HOME=${NGS_NXF_HOME:-"$HOME/data/.nextflow"}
export NXF_TEMP=${NGS_NXF_TEMP:-"/local/scratch/$USER/nextflow-preflight"}
export NXF_SYNTAX_PARSER=v2
mkdir -p "$NXF_HOME" "$NXF_TEMP"
cd "$repo_dir"
nextflow lint main.nf
nextflow run . \
  -profile hellbender_single_node \
  -params-file "$config_file" \
  -preview \
  -work-dir "$NXF_TEMP/work"
printf 'Nextflow preflight passed. No analysis tasks were executed.\n'
