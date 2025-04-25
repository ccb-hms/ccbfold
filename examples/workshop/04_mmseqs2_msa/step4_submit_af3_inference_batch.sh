#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=gpu_quad     # Partition to run in
#SBATCH -o logs/af3_inference_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/af3_inference_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM
#SBATCH --array=0-9              # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file
INPUT_DIR="af3_inference_inputs"
INPUT_FILES=(${INPUT_DIR}/*.json)
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

# Derive job name and output directory from file name
JOB_NAME=$(basename "$INPUT_FILE" .json)
OUTPUT_DIR="af3_inference_outputs"
JOB_OUTPUT_DIR="${OUTPUT_DIR}/${JOB_NAME}"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$JOB_OUTPUT_DIR"
