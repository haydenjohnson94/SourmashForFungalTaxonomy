#!/bin/bash
#SBATCH --job-name=DownloadGenomes_Cryptococcus
#SBATCH --output=DownloadGenomes_Cryptococcus_%j.out
#SBATCH --error=DownloadGenomes_Cryptococcus_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Downloads the Cryptococcus/Tremella species-complex genome assemblies used in
# this analysis. The accession list is extracted directly from the BUSCO tree
# shipped in this repo, rather than hardcoded, so it stays in sync with the
# data (leaf labels look like "Organism Name (GCA_XXXXXXXXX.X)").

DATASETS_TOOL="/path/to/your/datasets_tool/datasets"
TREE_FILE="../../SupplementalData_Accessions_Trees_And_BUSCOdatabases/CompareSourmashAndPhylogeneticAnalysis/Complexes/Cryptococcus/Cryptococcus_BUSCO.tree"
OUTPUT_DIR="/path/to/your/analysis_directory/SpeciesComplex_SourmashTrees/CryptococcusComplexes_Reduced/BUSCO_CryptococcusReduced"
DOWNLOAD_DIR="$OUTPUT_DIR/downloads"
GENOME_DIR="$OUTPUT_DIR/genomes"

mkdir -p "$DOWNLOAD_DIR" "$GENOME_DIR"

# Extract unique accessions (e.g. GCA_002215765.1) from the tree's leaf labels
mapfile -t ACCESSIONS < <(grep -oE "GCA_[0-9]+\.[0-9]+" "$TREE_FILE" | sort -u)
echo "Found ${#ACCESSIONS[@]} accessions in $TREE_FILE"

for ACCESSION in "${ACCESSIONS[@]}"; do
    echo "=== $ACCESSION ==="
    GENOME_COPY="$GENOME_DIR/${ACCESSION}_genomic.fna"

    if [[ -f "$GENOME_COPY" ]]; then
        echo "  Genome already present, skipping."
        continue
    fi

    ZIP_PATH="$DOWNLOAD_DIR/${ACCESSION}.zip"
    "$DATASETS_TOOL" download genome accession "$ACCESSION" --filename "$ZIP_PATH" --exclude-atypical

    if [[ ! -f "$ZIP_PATH" ]]; then
        echo "  ERROR: Download failed for $ACCESSION, skipping."
        continue
    fi

    unzip -o "$ZIP_PATH" -d "$DOWNLOAD_DIR/$ACCESSION"
    FNA_FILE=$(find "$DOWNLOAD_DIR/$ACCESSION" -name "*_genomic.fna" | head -1)

    if [[ -z "$FNA_FILE" ]]; then
        echo "  ERROR: No *_genomic.fna found for $ACCESSION after extraction, skipping."
        rm -rf "$DOWNLOAD_DIR/$ACCESSION" "$ZIP_PATH"
        continue
    fi

    mv "$FNA_FILE" "$GENOME_COPY"
    rm -rf "$DOWNLOAD_DIR/$ACCESSION" "$ZIP_PATH"
    echo "  Downloaded and extracted -> $GENOME_COPY"
done

echo "Download complete: $(ls "$GENOME_DIR" | wc -l) genomes in $GENOME_DIR"
