#!/bin/bash
#SBATCH --job-name=AF3Workshop_02_step2
#SBATCH -c 20                      # Request 20 cores
#SBATCH --mem=64G                  # Memory total in GiB
#SBATCH --partition=gpu            # Partition to run in
#SBATCH -o logs/af3_job_%A_%a.out  # Output log file, job ID (%A), task ID (%a)
#SBATCH --gres=gpu:l40s:1          # GPU requested
#SBATCH -t 0-01:00                 # Runtime in D-HH:MM
#SBATCH --array=0-9                # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file
INPUT_DIR="af3_inputs"
INPUT_FILES=(${INPUT_DIR}/*.json)
INPUT_FILE=${INPUT_FILES[$SLURM_ARRAY_TASK_ID]}

# Make sure output directory exists
OUTPUT_DIR="af3_outputs"
mkdir -p "$OUTPUT_DIR"

run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$OUTPUT_DIR"
