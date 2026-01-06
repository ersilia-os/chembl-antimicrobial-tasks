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
export CUDA_VISIBLE_DEVICES=0

# 🔍 GPU check *inside the container*
echo "=== Inside container GPU check ==="
echo "PATH=$PATH"
command -v nvidia-smi || echo "nvidia-smi NOT found in PATH"
nvidia-smi || echo "nvidia-smi FAILED inside container"

which ollama
echo "alpha=$alpha"
echo "OLLAMA_HOST=$OLLAMA_HOST"
echo "OLLAMA_MODELS=$OLLAMA_MODELS"
echo "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES"

# Start Ollama server and log directly to scratch (no cp at the end)
OLLAMA_LOG="/aloy/scratch/acomajuncosa/Ersilia/chembl_antimicrobial_tasks/ollama_logs/ollama_serve_${alpha}.log"
mkdir -p "$(dirname "$OLLAMA_LOG")"

echo "Starting ollama serve, logging to $OLLAMA_LOG"
ollama serve 2>&1 | tee "$OLLAMA_LOG" &
sleep 15

# Check whether it is responding
echo "Curling tags on $OLLAMA_HOST ..."
curl "http://127.0.0.1:${PORT}/api/tags" || echo "curl to ollama failed"

# Determine the directory where thisfile is located
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

# Call main.py using an absolute path
python "$SCRIPT_DIR/08_get_assay_descriptions_III.py" "$alpha"


