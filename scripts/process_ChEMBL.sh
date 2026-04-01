#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

#python $SCRIPT_DIR/00_export_chembl_activities.py
#python $SCRIPT_DIR/01_prepare_manual_files.py

echo ""
echo "ACTION REQUIRED before continuing:"
echo "  Fill in manual_curation_direction in:"
echo "    data/chembl_processed/01_activity_std_units_converted.csv"
echo "  Save the result as:"
echo "    config/activity_std_units_manual_curation.csv"
echo ""

CURATION_FILE="$SCRIPT_DIR/../config/activity_std_units_manual_curation.csv"
if [[ ! -f "$CURATION_FILE" ]]; then
  echo "ERROR: $CURATION_FILE not found. Complete manual curation before proceeding." >&2
  exit 1
fi

#python $SCRIPT_DIR/02_get_compound_info.py
python $SCRIPT_DIR/03_standardize_compounds.py
python $SCRIPT_DIR/04_merge_activity_and_compounds.py
python $SCRIPT_DIR/05_clean_activities.py


CALCULATE_ECFPS=0
for arg in "$@"; do
  [[ "$arg" == "--calculate_ecfps" ]] && CALCULATE_ECFPS=1 || {
    echo "Unknown argument: $arg" >&2
    exit 2
  }
done

if (( CALCULATE_ECFPS )); then
  python $SCRIPT_DIR/06_calculate_ecfps.py
fi