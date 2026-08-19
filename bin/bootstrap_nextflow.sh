#!/usr/bin/env bash
# Build the Nextflow controller and optional Shiny environments on a compute node.

set -euo pipefail
umask 0027

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/.." && pwd)
nextflow_prefix=${NGS_NEXTFLOW_ENV:-"$HOME/data/conda/envs/hellbender-ngs-nextflow"}
shiny_prefix=${NGS_SHINY_ENV:-"$HOME/data/conda/envs/hellbender-ngs-shiny"}
with_shiny=false

usage() {
  printf 'Usage: %s [--with-shiny] [--nextflow-prefix PATH] [--shiny-prefix PATH]\n' "$0"
}

while (($#)); do
  case "$1" in
    --with-shiny)
      with_shiny=true
      shift
      ;;
    --nextflow-prefix)
      nextflow_prefix=${2:?Missing value for --nextflow-prefix}
      shift 2
      ;;
    --shiny-prefix)
      shiny_prefix=${2:?Missing value for --shiny-prefix}
      shift 2
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

if [[ -z ${SLURM_JOB_ID:-} ]] && [[ $(hostname -s) == hellbender-* ]]; then
  printf 'ERROR: Run bootstrap inside an interactive compute allocation.\n' >&2
  printf 'Example: srun -p interactive --time=02:00:00 --mem=16G --cpus-per-task=4 --pty bash\n' >&2
  exit 2
fi

if command -v module >/dev/null 2>&1; then
  module purge
  module load miniconda3
fi
if ! command -v conda >/dev/null 2>&1; then
  printf 'ERROR: conda is unavailable. On Hellbender, run: module load miniconda3\n' >&2
  exit 2
fi

conda_base=$(conda info --base)
# The module determines this path at runtime.
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"

create_or_update() {
  local prefix=$1
  local environment_file=$2
  mkdir -p "$(dirname -- "$prefix")"
  if [[ -d $prefix/conda-meta ]]; then
    conda env update --prefix "$prefix" --file "$environment_file" --prune
  else
    conda env create --prefix "$prefix" --file "$environment_file"
  fi
}

printf 'Creating/updating Nextflow environment: %s\n' "$nextflow_prefix"
create_or_update "$nextflow_prefix" "$repo_dir/envs/nextflow.yaml"

if [[ $with_shiny == true ]]; then
  printf 'Creating/updating Shiny environment: %s\n' "$shiny_prefix"
  create_or_update "$shiny_prefix" "$repo_dir/envs/shiny.yaml"
fi

conda activate "$nextflow_prefix"
export NXF_SYNTAX_PARSER=v2
nextflow -version
printf 'Nextflow environment ready. Export this only for a non-default path:\n'
printf '  export NGS_NEXTFLOW_ENV=%q\n' "$nextflow_prefix"
if [[ $with_shiny == true ]]; then
  printf '  export NGS_SHINY_ENV=%q\n' "$shiny_prefix"
fi
