#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

if [[ -z "${1:-}" ]]; then
  echo "Usage: $0 --<pathogen_code>"
  exit 1
fi
code="${1#--}"

python "$SCRIPT_DIR/07_get_pathogen_assays.py" "$code"
python "$SCRIPT_DIR/08_clean_pathogen_activities.py" "$code"
python "$SCRIPT_DIR/09_curate_assay_parameters.py" "$code"
python "$SCRIPT_DIR/10_calculate_assay_clusters.py" "$code"
python "$SCRIPT_DIR/11_get_assay_overlap.py" "$code"
python "$SCRIPT_DIR/12_prepare_assay_datasets.py" "$code"
python "$SCRIPT_DIR/13_lightmodel_individual.py" "$code"
python "$SCRIPT_DIR/14_select_datasets_individual.py" "$code"
python "$SCRIPT_DIR/15_lightmodel_merged.py" "$code"
python "$SCRIPT_DIR/16_select_datasets_merged.py" "$code"
python "$SCRIPT_DIR/17_evaluate_correlations.py" "$code"
python "$SCRIPT_DIR/18_prepare_assay_master.py" "$code"
python "$SCRIPT_DIR/19_diagnosis.py" "$code"
