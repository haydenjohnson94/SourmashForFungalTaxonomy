#!/bin/bash
#SBATCH --job-name=RunBUSCO_Hypocreaceae_SameSpecies
#SBATCH --output=RunBUSCO_SameSpecies_%j.out
#SBATCH --error=RunBUSCO_SameSpecies_%j.err
#SBATCH --mem=120G
#SBATCH --cpus-per-task=32
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Runs the initial BUSCO pass on each genome downloaded by DownloadGenomes_SameSpecies.sh.
# Run RunBUSCO_AutoDetectReruns.sh afterward to retry any runs that failed or were
# interrupted.

BUSCO_DB="/path/to/your/busco_downloads/lineages/fungi_odb10"
OUTPUT_DIR="/path/to/your/analysis_directory/FirstManuscriptTreeComparisons/Kmers_SpeciesGenusFamily/Hypocreaceae/BUSCO_Hypocreaceae"
GENOME_DIR="$OUTPUT_DIR/genomes"
BUSCO_DIR="$OUTPUT_DIR/busco_results"

mkdir -p "$BUSCO_DIR"

for GENOME_COPY in "$GENOME_DIR"/*_genomic.fna; do
    [[ -f "$GENOME_COPY" ]] || continue
    ACCESSION=$(basename "$GENOME_COPY" "_genomic.fna")
    echo "=== $ACCESSION ==="

    SUMMARY_COUNT=$(find "$BUSCO_DIR/busco_${ACCESSION}" -name "short_summary.*.txt" 2>/dev/null | wc -l)
    if [[ "$SUMMARY_COUNT" -gt 0 ]]; then
        echo "  BUSCO already completed, skipping."
        continue
    fi

    echo "  Running BUSCO..."
    cd "$BUSCO_DIR"
    busco -i "$GENOME_COPY" \
          -o "busco_${ACCESSION}" \
          -l "$BUSCO_DB" \
          -m genome \
          -c 32 \
          --tar \
          --force

    if [[ $? -eq 0 ]]; then
        echo "  ✓ BUSCO completed successfully for $ACCESSION"
    else
        echo "  ✗ BUSCO failed for $ACCESSION (run RunBUSCO_AutoDetectReruns.sh afterward to retry)"
    fi
done

echo "Initial BUSCO pass complete!"
