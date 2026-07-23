#!/bin/bash
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account
#SBATCH --job-name=Zip_Hypocreaceae_32
#SBATCH --output=Zip_Hypocreaceae_%j.out
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G

# Change to the directory containing all the BUSCO results
cd /path/to/your/analysis_directory/FirstManuscriptTreeComparisons/Kmers_SpeciesGenusFamily/Hypocreaceae/BUSCO_Hypocreaceae/busco_results

# Loop through all busco_GCA_* directories
for busco_dir in busco_GCA_*/; do
    # Remove trailing slash
    busco_dir="${busco_dir%/}"
    
    # Extract the accession number (everything after busco_)
    accession="${busco_dir#busco_}"
    
    # Check if run_fungi_odb10 directory exists
    if [[ -d "$busco_dir/run_fungi_odb10" ]]; then
        echo "Zipping $busco_dir/run_fungi_odb10 as run_fungi_odb10_${accession}.tar.gz using 32 cores"
        
        # Create tar.gz using pigz for parallel compression
        tar -I "pigz -p 32" -cf "run_fungi_odb10_${accession}.tar.gz" -C "$busco_dir" "run_fungi_odb10"
        
        # Verify the archive was created
        if [[ $? -eq 0 ]]; then
            echo "✓ Successfully created run_fungi_odb10_${accession}.tar.gz"
        else
            echo "✗ Failed to create run_fungi_odb10_${accession}.tar.gz"
        fi
    else
        echo "⚠ run_fungi_odb10 not found in $busco_dir"
    fi
done

echo "All done! Created archives:"
ls -la run_fungi_odb10_*.tar.gz