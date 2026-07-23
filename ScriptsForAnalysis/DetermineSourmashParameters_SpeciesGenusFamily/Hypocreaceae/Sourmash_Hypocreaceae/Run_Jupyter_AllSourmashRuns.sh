#!/bin/bash
#SBATCH --job-name=run_notebook_Hypocreaceae_SameSpecies_SourmashSigAndANI
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --output=run_notebook_Hypocreaceae_SameSpecies_SourmashSigAndANI.%j.out
#SBATCH --error=run_notebook_Hypocreaceae_SameSpecies_SourmashSigAndANI.%j.err
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Input Jupyter Notebook file
NOTEBOOK_FILE="SameSpecies_Hypocreaceae_ComputeSourmashSigAndANI.ipynb"
EXECUTED_NOTEBOOK="${NOTEBOOK_FILE%.ipynb}_executed1.ipynb"
# Execute the notebook
echo "Executing $NOTEBOOK_FILE..."
jupyter nbconvert --to notebook --execute "$NOTEBOOK_FILE" --output "$EXECUTED_NOTEBOOK"
