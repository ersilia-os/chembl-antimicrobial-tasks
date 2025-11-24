#!/bin/bash

# Initialize conda
eval "$(conda shell.bash hook)"
conda activate /aloy/home/acomajuncosa/anaconda3/envs/camt

# Get alpha from first argument
alpha=$1

# === Portable Ollama environment (container paths!) ===
export PATH=/aloy/home/acomajuncosa/programs/ollama:$PATH
export OLLAMA_HOME=/aloy/home/acomajuncosa/programs/ollama/ollama_lib
export OLLAMA_MODELS=/aloy/home/acomajuncosa/programs/ollama/ollama_models
PORT=$((11435 + alpha))
export OLLAMA_HOST="http://127.0.0.1:${PORT}"

which ollama
echo "alpha=$alpha"
echo "OLLAMA_HOST=$OLLAMA_HOST"
echo "OLLAMA_MODELS=$OLLAMA_MODELS"

# Start Ollama server
ollama serve >/tmp/ollama_serve.log 2>&1 &
sleep 10

# Check whether it is responding
curl "http://127.0.0.1:${PORT}/api/tags"

# Determine the directory where thisfile is located
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Call main.py using an absolute path
python "$SCRIPT_DIR/08_get_assay_descriptions_III.py" "$alpha"

# Save the serve log permanently
cp /tmp/ollama_serve.log "/aloy/scratch/acomajuncosa/Ersilia/chembl_antimicrobial_tasks/ollama_logs/ollama_serve_${alpha}.log"

