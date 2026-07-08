#!/usr/bin/env bash
# Generic helper: run one per-pathogen step across pathogens. Run in the `camt` env.
#
# Usage: run_step.sh <step_script> [pathogen_codes...]
#   e.g.: run_step.sh 10_calculate_assay_clusters.py
#         run_step.sh 21_build_curated_datasets.py mtuberculosis ecoli
#   default codes: all in config/pathogens.csv
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

if [[ -z "${1:-}" ]]; then
  echo "Usage: $0 <step_script> [pathogen_codes...]"
  echo "  e.g.: $0 10_calculate_assay_clusters.py"
  exit 1
fi

step="$1"; shift

mapfile -t CODES < <(tail -n +2 "$SCRIPT_DIR/../config/pathogens.csv" | cut -d, -f2)
[ "$#" -gt 0 ] && CODES=("$@")

for code in "${CODES[@]}"; do
  code="${code#--}"
  echo "=== $step $code ==="
  python "$SCRIPT_DIR/$step" "$code"
done
