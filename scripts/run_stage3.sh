#!/usr/bin/env bash
# Stage 3 - conservative individual-dataset branch (13-20), per pathogen. Run in `camt`.
# Parallel to Stage 4 (not a prerequisite for it); needs step 09 output (assay parameters).
# Opt-in: not part of run_all.sh.
#
# Usage: run_stage3.sh [pathogen_codes...]   (default: all codes in config/pathogens.csv)
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

mapfile -t CODES < <(tail -n +2 "$SCRIPT_DIR/../config/pathogens.csv" | cut -d, -f2)
[ "$#" -gt 0 ] && CODES=("$@")

STEPS=(
  13_prepare_assay_datasets
  14_lightmodel_individual
  15_select_datasets_individual
  16_lightmodel_merged
  17_select_datasets_merged
  18_evaluate_correlations
  19_prepare_assay_master
  20_diagnosis
)

for code in "${CODES[@]}"; do
  code="${code#--}"
  echo "=== Stage 3 - $code ==="
  for step in "${STEPS[@]}"; do
    python "$SCRIPT_DIR/${step}.py" "$code"
  done
done
