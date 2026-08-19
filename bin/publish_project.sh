#!/usr/bin/env bash
# Publish a reviewed project through a symlink without changing permissions.

set -euo pipefail
umask 0027

if (($# != 2)); then
  printf 'Usage: %s COMPLETED_PROJECT_DIR PORTAL_ROOT\n' "$0" >&2
  exit 2
fi

project_dir=$1
portal_root=$2
if [[ ! -f $project_dir/portal/manifest.json ]]; then
  printf 'ERROR: Completed portal manifest not found: %s\n' "$project_dir/portal/manifest.json" >&2
  exit 2
fi

project_dir=$(cd -- "$project_dir" && pwd)
mkdir -p "$portal_root"
portal_root=$(cd -- "$portal_root" && pwd)

project_id=$(python3 - "$project_dir/portal/manifest.json" <<'PY'
import json
import re
import sys

with open(sys.argv[1], encoding="utf-8") as handle:
    project_id = json.load(handle)["project"]["id"]
if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9._-]{0,99}", project_id):
    raise SystemExit("Unsafe project id in manifest")
print(project_id)
PY
)

destination=$portal_root/$project_id
if [[ -L $destination ]]; then
  current_target=$(readlink -f -- "$destination")
  if [[ $current_target == "$project_dir" ]]; then
    printf 'Project is already published: %s\n' "$destination"
    exit 0
  fi
  printf 'ERROR: Destination already points to another project: %s\n' "$destination" >&2
  exit 2
fi
if [[ -e $destination ]]; then
  printf 'ERROR: Destination already exists: %s\n' "$destination" >&2
  exit 2
fi

ln -s -- "$project_dir" "$destination"
printf 'Published %s at %s\n' "$project_id" "$destination"
printf 'No permissions were changed. Verify the intended client group can read the target.\n'

