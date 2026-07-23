#!/bin/bash
#SBATCH --job-name=run_notebook_BuildBUSCOTree_Cryptococcus
#SBATCH --cpus-per-task=32
#SBATCH --mem=150G
#SBATCH --output=run_notebook_BuildBUSCOTree_Cryptococcus.%j.out
#SBATCH --error=run_notebook_BuildBUSCOTree_Cryptococcus.%j.err
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Input Jupyter Notebook file
NOTEBOOK_FILE="BuildBUSCOTree_Cryptococcus_Reduced.ipynb"
EXECUTED_NOTEBOOK="${NOTEBOOK_FILE%.ipynb}_executed.ipynb"

# Execute the notebook
echo "Executing $NOTEBOOK_FILE..."
jupyter nbconvert --to notebook --execute "$NOTEBOOK_FILE" --output "$EXECUTED_NOTEBOOK"
