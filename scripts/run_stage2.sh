#!/usr/bin/env bash
# Stage 2 - per-pathogen assay preparation. Run in the `camt` env.
#
# Runs the camt-runnable per-pathogen steps: 07, 08, 10, 11. Two Stage-2 steps are NOT run
# here because they need a different environment - run them manually (see NOTICE below):
#   09_curate_assay_parameters.py  GPU + local ollama   (per pathogen; needed before Stages 3/4)
#   12_calculate_eos1klk.py        ersilia env          (cross-pathogen; needed for Stage-4 embeddings)
#
# Usage: run_stage2.sh [pathogen_codes...]   (default: all codes in config/pathogens.csv)
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

mapfile -t CODES < <(tail -n +2 "$SCRIPT_DIR/../config/pathogens.csv" | cut -d, -f2)
[ "$#" -gt 0 ] && CODES=("$@")

for code in "${CODES[@]}"; do
  code="${code#--}"
  echo "=== Stage 2 - $code ==="
  python "$SCRIPT_DIR/07_get_pathogen_assays.py" "$code"
  python "$SCRIPT_DIR/08_clean_pathogen_activities.py" "$code"
  python "$SCRIPT_DIR/10_calculate_assay_clusters.py" "$code"
  python "$SCRIPT_DIR/11_get_assay_overlap.py" "$code"
done

cat >&2 <<'NOTE'

------------------------------------------------------------------
Stage 2 NOTICE - two steps are NOT run here (they need other envs):
  09_curate_assay_parameters.py  (GPU + local ollama)  per pathogen
       REQUIRED before Stages 3 and 4
  12_calculate_eos1klk.py        (ersilia env)         cross-pathogen, run once
       REQUIRED for the Stage-4 embedding figures
Run them manually, e.g.:
  python scripts/09_curate_assay_parameters.py <code>          # ollama running
  conda run -n ersilia python scripts/12_calculate_eos1klk.py
------------------------------------------------------------------
NOTE
