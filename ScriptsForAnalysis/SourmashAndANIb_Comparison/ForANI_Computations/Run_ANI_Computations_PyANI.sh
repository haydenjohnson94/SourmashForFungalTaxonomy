#!/bin/bash
#SBATCH --job-name=ANI_Computations_PyANI
#SBATCH --output=ANI_Computations_PyANI_%j.out
#SBATCH --error=ANI_Computations_PyANI_%j.err
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Runs the PyANI ANIb pairwise-comparison pipeline for each evolutionary-path row.
# Requires: NCBI datasets tool, PyANI (average_nucleotide_identity.py) on PATH.
python ANI_Computations_PyANI.py
