#!/bin/bash

#SBATCH -J descriptions
#SBATCH --chdir=/aloy/home/acomajuncosa/Ersilia/chembl-antimicrobial-tasks
#SBATCH --ntasks=1
#SBATCH --array=0-1%5
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=/aloy/scratch/acomajuncosa/Ersilia/chembl_antimicrobial_tasks/descriptions/%x_%a.out
#SBATCH --gpus=1   
#SBATCH --time=300:00:00
#SBATCH -p sbnb_gpu_3090
#SBATCH -w irbgcn07  

# Loads default environment configuration
# export SINGULARITYENV_LD_LIBRARY_PATH=$LD_LIBRARY_PATH #/.singularity.d/libs
export SINGULARITY_BINDPATH="/home/sbnb:/aloy/home,/data/sbnb/data:/aloy/data,/data/sbnb/scratch:/aloy/scratch"

# Load cuda libraries
export LD_LIBRARY_PATH=/apps/manual/software/CUDA/11.6.1/lib64:/apps/manual/software/CUDA/11.6.1/targets/x86_64-linux/lib:/apps/manual/software/CUDA/11.6.1/extras/CUPTI/lib64/:/apps/manual/software/CUDA/11.6.1/nvvm/lib64/:$LD_LIBRARY_PATH

# Run jobs
alpha=($(seq 0 1))
singularity exec --cleanenv --nv /apps/singularity/ood_images/jupyter_ollama_0.11.5.sif ./scripts/08_get_assay_descriptions/08_get_assay_descriptions_II.sh ${alpha[$SLURM_ARRAY_TASK_ID]}