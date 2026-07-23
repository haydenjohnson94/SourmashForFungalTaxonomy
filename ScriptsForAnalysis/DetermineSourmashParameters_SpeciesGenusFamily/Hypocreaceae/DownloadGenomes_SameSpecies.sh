#!/bin/bash
#SBATCH --job-name=DownloadGenomes_Hypocreaceae_SameSpecies
#SBATCH --output=DownloadGenomes_SameSpecies_%j.out
#SBATCH --error=DownloadGenomes_SameSpecies_%j.err
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Downloads the same-species genome assemblies (all Trichoderma reesei; matches
# SymlinkAssemblies.sh) from NCBI into GENOME_DIR. Run RunBUSCO_SameSpecies.sh
# afterward to compute BUSCOs on the downloaded genomes.

DATASETS_TOOL="/path/to/your/datasets_tool/datasets"
OUTPUT_DIR="/path/to/your/analysis_directory/FirstManuscriptTreeComparisons/Kmers_SpeciesGenusFamily/Hypocreaceae/BUSCO_Hypocreaceae"
DOWNLOAD_DIR="$OUTPUT_DIR/downloads"
GENOME_DIR="$OUTPUT_DIR/genomes"

mkdir -p "$DOWNLOAD_DIR" "$GENOME_DIR"

# Same-species accessions (all Trichoderma reesei; matches SymlinkAssemblies.sh)
ACCESSIONS=("GCA_000167675.2" "GCA_002006585.1"
            "GCA_052818895.1" "GCA_047716485.1"
            "GCA_047716445.1" "GCA_047716435.1"
            "GCA_047716415.1" "GCA_042257965.1"
            "GCA_038428105.1" "GCA_028871515.1"
            "GCA_016806875.1" "GCA_028871095.1"
            "GCA_016806815.1" "GCA_028870775.1")

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
