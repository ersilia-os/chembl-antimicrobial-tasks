#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

# python $SCRIPT_DIR/../src/00_export_chembl_activities.py
# python $SCRIPT_DIR/../src/01_get_compound_info.py
# python $SCRIPT_DIR/../src/02_standardize_compounds.py
# python $SCRIPT_DIR/../src/03_merge_all.py
# python $SCRIPT_DIR/../src/04_prepare_conversions.py
# python $SCRIPT_DIR/../src/05_clean_activities.py


CALCULATE_ECFPS=0
for arg in "$@"; do
  [[ "$arg" == "--calculate_ecfps" ]] && CALCULATE_ECFPS=1 || {
    echo "Unknown argument: $arg" >&2
    exit 2
  }
done

if (( CALCULATE_ECFPS )); then
  python $SCRIPT_DIR/../src/06_calculate_ecfps.py
fi