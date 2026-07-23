#!/bin/bash
#SBATCH --job-name=ANI_Computations_Sourmash
#SBATCH --output=ANI_Computations_Sourmash_%j.out
#SBATCH --error=ANI_Computations_Sourmash_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Runs the sourmash signature/ANI pairwise-comparison pipeline for each evolutionary-path row.
# Requires: NCBI datasets tool, sourmash on PATH.
python ANI_Computations_Sourmash.py
