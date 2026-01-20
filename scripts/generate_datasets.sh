#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

if [[ -z "${1:-}" ]]; then
  echo "Usage: $0 --<pathogen_code>"
  exit 1
fi
code="${1#--}"

python "$SCRIPT_DIR/../src/07_get_pathogen_assays.py" "$code"
# python "$SCRIPT_DIR/../src/08_clean_pathogen_activities.py" "$code"
