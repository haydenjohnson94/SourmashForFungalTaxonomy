#!/bin/bash
#SBATCH --job-name=BUSCO_Rerun_Hypocreaceae
#SBATCH --output=BUSCO_Rerun_%j.out
#SBATCH --error=BUSCO_Rerun_%j.err
#SBATCH --mem=120G
#SBATCH --cpus-per-task=16
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Run this after DownloadGenomes_SameSpecies.sh and RunBUSCO_SameSpecies.sh, to
# retry any BUSCO runs that failed or were interrupted during the initial pass.

# Paths
BUSCO_DB="/path/to/your/busco_downloads/lineages/fungi_odb10"
OUTPUT_DIR="/path/to/your/analysis_directory/FirstManuscriptTreeComparisons/Kmers_SpeciesGenusFamily/Hypocreaceae/BUSCO_Hypocreaceae"
GENOME_DIR="$OUTPUT_DIR/genomes"
BUSCO_DIR="$OUTPUT_DIR/busco_results"

# --- DETECTION LOGIC ---
echo "Scanning for failed BUSCO runs in $BUSCO_DIR..."
RERUN_ACCESSIONS=()

# Check every directory starting with "busco_GCA_"
for dir in "$BUSCO_DIR"/busco_GCA_*; do
    if [ -d "$dir" ]; then
        ACCESSION=$(basename "$dir" | sed 's/busco_//')
        
        # Check if directory is empty OR missing the short_summary file
        # (This is more robust than just checking if it's empty)
        SUMMARY_COUNT=$(find "$dir" -name "short_summary.*.txt" | wc -l)
        
        if [ ! "$(ls -A "$dir")" ] || [ "$SUMMARY_COUNT" -eq 0 ]; then
            echo "  Detected failure: $ACCESSION (Empty or no summary)"
            RERUN_ACCESSIONS+=("$ACCESSION")
        fi
    fi
done

if [ ${#RERUN_ACCESSIONS[@]} -eq 0 ]; then
    echo "No failed runs detected. Exiting."
    exit 0
fi

echo "Found ${#RERUN_ACCESSIONS[@]} accessions to rerun."
echo "-----------------------------------------------"

# --- RERUN LOOP ---
for ACCESSION in "${RERUN_ACCESSIONS[@]}"; do
    echo "Processing Rerun for $ACCESSION..."
    
    # Path to the genome we previously extracted
    GENOME_COPY="$GENOME_DIR/${ACCESSION}_genomic.fna"
    
    if [[ -f "$GENOME_COPY" ]]; then
        cd "$BUSCO_DIR"
        echo "  Running BUSCO on $ACCESSION with --force..."
        
        # Using --force (or -f) to overwrite the existing folder
        busco -i "$GENOME_COPY" \
              -o "busco_${ACCESSION}" \
              -l "$BUSCO_DB" \
              -m genome \
              -c 32 \
              --tar \
              --force
        
        if [[ $? -eq 0 ]]; then
            echo "  ✓ BUSCO rerun completed successfully for $ACCESSION"
        else
            echo "  ✗ BUSCO rerun failed again for $ACCESSION"
        fi
    else
        # If the genome file isn't in GENOME_DIR, re-run DownloadGenomes_SameSpecies.sh
        echo "  ✗ Genome file $GENOME_COPY not found. Skipping."
    fi
    echo ""
done

echo "Rerun analysis complete!"