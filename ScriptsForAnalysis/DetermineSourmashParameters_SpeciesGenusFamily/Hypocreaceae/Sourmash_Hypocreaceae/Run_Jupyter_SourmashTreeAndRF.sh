#!/bin/bash
#SBATCH --job-name=run_notebook_Hypocreaceae_SameSpecies_SourmashTreeAndRF
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --output=run_notebook_Hypocreaceae_SameSpecies_SourmashTreeAndRF.%j.out
#SBATCH --error=run_notebook_Hypocreaceae_SameSpecies_SourmashTreeAndRF.%j.err
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Run this after Run_Jupyter_AllSourmashRuns.sh (which computes the sourmash
# signatures/ANI matrices this notebook builds trees from and computes RF
# distance against a reference tree).

# Input Jupyter Notebook file
NOTEBOOK_FILE="SameSpecies_Hypocreaceae_SourmashTreeComparisons_UPGMAtrees_And_RF_Computation.ipynb"
EXECUTED_NOTEBOOK="${NOTEBOOK_FILE%.ipynb}_executed.ipynb"

# Execute the notebook
echo "Executing $NOTEBOOK_FILE..."
jupyter nbconvert --to notebook --execute "$NOTEBOOK_FILE" --output "$EXECUTED_NOTEBOOK"
