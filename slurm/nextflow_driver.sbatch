#!/usr/bin/env bash
#SBATCH --job-name=ngs-nextflow
#SBATCH --partition=general
#SBATCH --time=2-00:00:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --export=ALL

set -euo pipefail
umask 0027

: "${NGS_CONFIG:?NGS_CONFIG was not exported by submit_nextflow.sh}"
: "${NGS_NXF_PROFILE:?NGS_NXF_PROFILE was not exported by submit_nextflow.sh}"
: "${NGS_NEXTFLOW_ENV:?NGS_NEXTFLOW_ENV was not exported by submit_nextflow.sh}"
: "${NGS_REPO_DIR:?NGS_REPO_DIR was not exported by submit_nextflow.sh}"

module purge
module load miniconda3
conda_base=$(conda info --base)
# The module determines this path at runtime.
# shellcheck disable=SC1091
source "$conda_base/etc/profile.d/conda.sh"
conda activate "$NGS_NEXTFLOW_ENV"

config_values=$(python3 - "$NGS_CONFIG" <<'PY'
import sys
from pathlib import Path
import yaml

with Path(sys.argv[1]).open(encoding="utf-8") as handle:
    config = yaml.safe_load(handle) or {}
output = Path(config["project"]["output_dir"]).expanduser().resolve()
project_id = config["project"]["id"]
print(output)
print(project_id)
PY
)
output_dir=$(sed -n '1p' <<<"$config_values")
project_id=$(sed -n '2p' <<<"$config_values")
work_dir=$output_dir/work/nextflow
provenance_dir=$output_dir/provenance
mkdir -p "$work_dir" "$provenance_dir"
attempt_id=${SLURM_JOB_ID}.${SLURM_RESTART_COUNT:-0}

export NXF_HOME=${NGS_NXF_HOME:-"$HOME/data/.nextflow"}
export NXF_SYNTAX_PARSER=v2
driver_tmp=/local/scratch/$USER/nextflow-driver-${SLURM_JOB_ID}
mkdir -p "$NXF_HOME" "$driver_tmp"
export NXF_TEMP=$driver_tmp
trap 'rm -rf -- "$driver_tmp"' EXIT

cd "$NGS_REPO_DIR"
python workflow/lib/pipeline_config.py --config "$NGS_CONFIG"
printf 'Controller job: %s\n' "$SLURM_JOB_ID"
nextflow_version=$(nextflow -version | awk '
  /version/ {
    for (i = 1; i <= NF; i++) {
      if ($i == "version") {
        print $(i + 1)
        exit
      }
    }
  }
')
printf 'Nextflow: %s\n' "$nextflow_version"
printf 'Project: %s\n' "$project_id"
printf 'Output: %s\n' "$output_dir"

nextflow_args=(
  run "$NGS_REPO_DIR"
  -profile "$NGS_NXF_PROFILE"
  -params-file "$NGS_CONFIG"
  -work-dir "$work_dir"
  -resume
  -with-trace "$provenance_dir/nextflow_trace.${attempt_id}.tsv"
  -with-report "$provenance_dir/nextflow_report.${attempt_id}.html"
  -with-timeline "$provenance_dir/nextflow_timeline.${attempt_id}.html"
  -with-dag "$provenance_dir/nextflow_dag.${attempt_id}.html"
)
nextflow_command=(nextflow)
if [[ -n ${NGS_NXF_CONFIG:-} ]]; then
  # -c is a Nextflow launcher option and must precede the `run` command.
  nextflow_command+=(-c "$NGS_NXF_CONFIG")
fi
nextflow_command+=("${nextflow_args[@]}")
"${nextflow_command[@]}"
