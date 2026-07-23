#!/bin/bash
#SBATCH --job-name=RunBUSCO_Cryptococcus
#SBATCH --output=RunBUSCO_Cryptococcus_%j.out
#SBATCH --error=RunBUSCO_Cryptococcus_%j.err
#SBATCH --mem=900G
#SBATCH --cpus-per-task=32
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account

# Runs the initial BUSCO pass on each genome downloaded by DownloadGenomes_Cryptococcus.sh.
# Afterward, run Zip_BUSCO_Results_Cryptococcus.sh to package results in the layout
# BuildBUSCOTree_Cryptococcus_Reduced.ipynb expects.

BUSCO_DB="/path/to/your/busco_downloads/lineages/fungi_odb10"
OUTPUT_DIR="/path/to/your/analysis_directory/SpeciesComplex_SourmashTrees/CryptococcusComplexes_Reduced/BUSCO_CryptococcusReduced"
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
        echo "  ✗ BUSCO failed for $ACCESSION"
    fi
done

echo "Initial BUSCO pass complete!"
