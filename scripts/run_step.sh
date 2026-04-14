#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)

if [[ -z "${1:-}" ]]; then
  echo "Usage: $0 <step_script>"
  echo "  e.g.: $0 10_calculate_assay_clusters.py"
  exit 1
fi

step="$1"

pathogens=(
  abaumannii
  mtuberculosis
  campylobacter
  ngonorrhoeae
  hpylori
  smansoni
  enterobacter
  efaecium
  spneumoniae
  kpneumoniae
  pfalciparum
  calbicans
  paeruginosa
  ecoli
  saureus
)

for code in "${pathogens[@]}"; do
  echo "=== $step $code ==="
  python "$SCRIPT_DIR/$step" "$code"
done
