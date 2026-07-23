#!/bin/bash
#SBATCH --job-name=DataframeBuilding_Analysis
#SBATCH --output=DataframeBuildingOutput_%j.out
#SBATCH --error=DataframeBuildingError_%j.err
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Run the accessions-dataframe script
python /path/to/your/analysis_directory/KingdomToSpeciesANI/SpeciesGenusAnalysis/DetermineAccessionsDataframe_SpeciesGenus.py
