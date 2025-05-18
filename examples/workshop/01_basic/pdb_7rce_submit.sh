#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=gpu          # Partition to run in
#SBATCH -o pdb_7rce_%j.out       # Output log file, job ID (%j)
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM format

module load alphafold/3.0.1

run_alphafold.py \
   --json_path=pdb_7rce_input.json \
   --output_dir=.