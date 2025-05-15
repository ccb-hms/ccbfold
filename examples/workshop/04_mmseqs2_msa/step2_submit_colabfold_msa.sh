#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=gpu          # Partition to run in
#SBATCH -o logs/colabfold_msa_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/colabfold_msa_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-08:00               # Runtime in D-HH:MM

module load localcolabfold/406d4c6
module load mmseqs2/17-b804f

OUTPUT_DIR="colabfold_msa_outputs"

# Shared database directory for colabfold MSAs
DATABASE_DIR="/n/shared_db/colabfold/406d4c6/"

# Run MSAs with mmseqs2 on GPU
colabfold_search \
   --gpu 1 \
   ../data/preliminaryTests.fasta \
   $DATABASE_DIR \
   $OUTPUT_DIR