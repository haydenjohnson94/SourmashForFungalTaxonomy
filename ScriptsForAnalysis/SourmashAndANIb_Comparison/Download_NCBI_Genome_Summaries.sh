#!/bin/bash

# Define the tool path
DATASETS="/path/to/your/datasets_tool/datasets"

# Date stamp used to version these downloads (matches the naming convention expected
# by DetermineEvolutionaryPaths.py and DetermineAccessionsDataframe_SpeciesGenus.py,
# e.g. fungi_genome_summary_Apr14_2025.json). Update those scripts' filenames to match
# if you regenerate this on a new date.
DATE_STAMP=$(date +%b%d_%Y)

# NCBI's datasets tool outputs JSON by default
$DATASETS summary genome taxon 'fungi' > "fungi_genome_summary_${DATE_STAMP}.json"
$DATASETS summary genome taxon 'oomycota' > "oomycete_genome_summary_${DATE_STAMP}.json"