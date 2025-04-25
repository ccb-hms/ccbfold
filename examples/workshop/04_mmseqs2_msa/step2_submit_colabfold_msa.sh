#!/bin/bash
#SBATCH -c 20                    # Request 20 cores
#SBATCH --mem=64G                # Memory total in GiB
#SBATCH --partition=short        # Partition to run in
#SBATCH -o logs/colabfold_msa_job_%A_%a.out       # STDOUT file
#SBATCH -e logs/colabfold_msa_job_%A_%a.err       # STDERR file
#SBATCH --gres=gpu:l40s:1        # GPU requested
#SBATCH -t 0-01:00               # Runtime in D-HH:MM

module use /n/shared_db/tmp/module
module load colabfold
module load mmseqs2

OUTPUT_DIR="colabfold_msa_outputs"

# Shared database directory for colabfold MSA
DATABASE_DIR="/n/shared_db/tmp/colabfold"

# Run MSA with mmseqs2 on GPU
colabfold_search \
   --gpu 1 \
   ../data/preliminaryTests.fasta \
   $DATABASE_DIR \
   $OUTPUT_DIR