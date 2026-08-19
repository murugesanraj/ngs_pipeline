#!/usr/bin/env bash
# Create private project configuration from tracked examples.

set -euo pipefail
umask 0027

if (($# != 1)); then
  printf 'Usage: %s PROJECT_DIRECTORY\n' "$0" >&2
  exit 2
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
repo_dir=$(cd -- "$script_dir/.." && pwd)
project_dir=$1

if [[ $project_dir == *$'\n'* || $project_dir == *,* ]]; then
  printf 'ERROR: Project path cannot contain a newline or comma.\n' >&2
  exit 2
fi

mkdir -p "$project_dir"
project_dir=$(cd -- "$project_dir" && pwd)

for filename in config.yaml samples.tsv contrasts.tsv; do
  if [[ -e $project_dir/$filename ]]; then
    printf 'ERROR: Refusing to overwrite %s\n' "$project_dir/$filename" >&2
    exit 2
  fi
done

cp "$repo_dir/config/config.example.yaml" "$project_dir/config.yaml"
cp "$repo_dir/config/samples.example.tsv" "$project_dir/samples.tsv"
cp "$repo_dir/config/contrasts.example.tsv" "$project_dir/contrasts.tsv"

python3 - "$project_dir/config.yaml" "$project_dir" <<'PY'
from pathlib import Path
import sys

config_path = Path(sys.argv[1])
project_path = sys.argv[2]
text = config_path.read_text(encoding="utf-8")
text = text.replace(
    "/cluster/pixstor/GROUP/ngs-projects/example_project",
    project_path,
)
config_path.write_text(text, encoding="utf-8")
PY

printf 'Project configuration created in %s\n' "$project_dir"
printf 'Edit these files before validation:\n'
printf '  %s\n' "$project_dir/config.yaml" "$project_dir/samples.tsv" "$project_dir/contrasts.tsv"

