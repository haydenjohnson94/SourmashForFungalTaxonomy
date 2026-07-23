This repository contains supplemental data for the manuscript "Whole genome similarity provides a robust framework for classification of fungal taxa from the genus rank to infraspecific variants "

## Repository layout

- **`SupplementalData_Accessions_Trees_And_BUSCOdatabases/`** - tree files, genome
  accession lists, and other supporting data from the manuscript itself.

- **`ScriptsForAnalysis/`** - scripts and Jupyter notebooks to recreate the types of
  analysis done in the manuscript (BUSCO- and sourmash-based tree building,
  Robinson-Foulds tree comparisons, sourmash-vs-ANIb evolutionary path comparisons).
  Each subdirectory has its own `RunOrder.txt` explaining the order to run its
  scripts in.

- **`Tutorial_QuickStart_ANItree/`** - a short, self-contained walkthrough of the
  core workflow (download genomes &rarr; sketch sourmash signatures &rarr; compute
  pairwise ANI &rarr; heatmap &rarr; UPGMA tree) on a single example genus, with real
  executed outputs included in the notebook. Start here if you just want to see the
  method work end-to-end before diving into the full analyses above.

## Data availability note

The NCBI genome-catalog JSON snapshot used as input to the evolutionary-path analysis
(`ScriptsForAnalysis/SourmashAndANIb_Comparison/SpeciesGenusAnalysis/fungi_genome_summary_Apr14_2025.json`)
is not reproducibly regeneratable (NCBI's catalog grows over time) and is excluded
from this repository via `.gitignore`. It is archived on Zenodo instead
(DOI: TODO - add link here once the Zenodo deposit is live).
