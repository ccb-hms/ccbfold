#!/bin/bash
#SBATCH --job-name=AF3Workshop_03_step3
#SBATCH -c 2                                 # Request 2 cores
#SBATCH --mem=16G                            # Memory total in GiB
#SBATCH --partition=gpu                      # Partition to run in
#SBATCH -o logs/af3_inference_job_%A_%a.out  # Output log file
#SBATCH --gres=gpu:l40s:1                    # GPU requested
#SBATCH -t 0-01:00                           # Runtime in D-HH:MM
#SBATCH --array=0-9                          # Job array indices (for 10 .json files)

module load alphafold/3.0.1

# Get the input json file with MSA created in step 2
INPUT_DIR="af3_msa_outputs"
JOB_NAMES=($(ls "${INPUT_DIR}"))
JOB_NAME=${JOB_NAMES[$SLURM_ARRAY_TASK_ID]}
INPUT_FILE="${INPUT_DIR}/${JOB_NAME}/${JOB_NAME}_data.json"

# Make sure output directory exists
OUTPUT_DIR="af3_inference_outputs"
mkdir -p "$OUTPUT_DIR"

# Run AlphaFold3 Inference only
run_alphafold.py \
   --json_path="$INPUT_FILE" \
   --output_dir="$OUTPUT_DIR" \
   --norun_data_pipeline
