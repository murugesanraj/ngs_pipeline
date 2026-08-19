#!/usr/bin/env bash
#SBATCH --job-name=ngs-driver
#SBATCH --partition=general
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL

set -euo pipefail
umask 0027

: "${NGS_CONFIG:?NGS_CONFIG was not exported by submit_pipeline.sh}"
: "${NGS_PROFILE:?NGS_PROFILE was not exported by submit_pipeline.sh}"
: "${NGS_WORKFLOW_ENV:?NGS_WORKFLOW_ENV was not exported by submit_pipeline.sh}"
: "${NGS_REPO_DIR:?NGS_REPO_DIR was not exported by submit_pipeline.sh}"

module purge
module load miniconda3
conda_base=$(conda info --base)
# The module determines this path at runtime.
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
conda activate "$NGS_WORKFLOW_ENV"

cd "$NGS_REPO_DIR"
printf 'Controller job: %s\n' "${SLURM_JOB_ID:-unknown}"
printf 'Host: %s\n' "$(hostname -f)"
printf 'Snakemake: %s\n' "$(snakemake --version)"
printf 'Config: %s\n' "$NGS_CONFIG"
printf 'Profile: %s\n' "$NGS_PROFILE"

snakemake \
  --snakefile "$NGS_REPO_DIR/workflow/Snakefile" \
  --configfile "$NGS_CONFIG" \
  --profile "$NGS_PROFILE"
