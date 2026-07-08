#!/usr/bin/env bash
# Full repeatable run: Stage 2 -> Stage 4, over all pathogens in config/pathogens.csv
# (pass codes to restrict, e.g. `run_all.sh mtuberculosis ecoli`). Run in the `camt` env.
#
# PREREQUISITES (not run here):
#   - Stage 1, once:   bash run_stage1.sh   [--calculate_ecfps]
#   - Step 09 (GPU + ollama, per pathogen) and step 12 (ersilia, once) - see run_stage2.sh NOTICE.
#     09 is required before Stage 4; 12 is required for the Stage-4 embedding figures.
#   - Stage 3 (conservative branch) is a separate opt-in run:  bash run_stage3.sh
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

bash "$SCRIPT_DIR/run_stage2.sh" "$@"
bash "$SCRIPT_DIR/run_stage4.sh" "$@"
