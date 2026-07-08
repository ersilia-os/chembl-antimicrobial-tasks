#!/usr/bin/env bash
# Stage 4 - signal-based pooling branch (21-28), per pathogen + general. Run in `camt`.
#
# Per-pathogen data block (dependency chain): 21 -> 22 -> 23 -> 24 -> 25.
# Then the general/reporting block once, over ALL pathogens on disk: 26 (low-data catch-all,
# loops pathogens internally), 27 (ALL cross-pathogen figures), 28 (dataset trace/flow audit).
# Embedding figures need step 12 output (skipped gracefully if absent).
#
# Usage: run_stage4.sh [pathogen_codes...]   (default: all codes in config/pathogens.csv)
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

mapfile -t CODES < <(tail -n +2 "$SCRIPT_DIR/../config/pathogens.csv" | cut -d, -f2)
[ "$#" -gt 0 ] && CODES=("$@")

for code in "${CODES[@]}"; do
  code="${code#--}"
  echo "=== Stage 4 (per-pathogen) - $code ==="
  python "$SCRIPT_DIR/21_build_curated_datasets.py" "$code"
  python "$SCRIPT_DIR/22_binarise_datasets.py" "$code"
  python "$SCRIPT_DIR/23_transfer_pooling.py" "$code"
  python "$SCRIPT_DIR/24_crossvalidate_pools.py" "$code"
  python "$SCRIPT_DIR/25_merge_leftovers.py" "$code"
done

# General/reporting block (once; always aggregates every pathogen on disk, independent of
# any code subset above - cross-pathogen figures want all available data).
echo "=== Stage 4 (general/reporting) ==="
python "$SCRIPT_DIR/26_lowdata_pools.py"
python "$SCRIPT_DIR/27_general_plots.py"
python "$SCRIPT_DIR/28_trace_sample_datasets.py"
