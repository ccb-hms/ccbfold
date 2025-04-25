#!/bin/bash
#SBATCH -c 10                    # Request 10 cores
#SBATCH --mem=16G                # Memory total in GiB
#SBATCH --partition=gpu_quad     # Partition to run in
#SBATCH -o logs/af3_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/af3_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM
#SBATCH --array=0-9              # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file with MSA created in step 2
INPUT_DIR="af3_msa_outputs"
JOB_NAMES=($(ls "${INPUT_DIR}"))
JOB_NAME=${JOB_NAMES[$SLURM_ARRAY_TASK_ID]}
INPUT_FILE="${INPUT_DIR}/${JOB_NAME}/${JOB_NAME}_data.json"

# Derive job name and output directory from file name
OUTPUT_DIR="af3_inference_outputs"
JOB_OUTPUT_DIR="${OUTPUT_DIR}/${JOB_NAME}"

# Make sure output directory exists
mkdir -p "$OUTPUT_DIR"

# Run AlphaFold3 Inference only
run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$JOB_OUTPUT_DIR" \
   --norun_data_pipeline
