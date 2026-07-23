#!/bin/bash
#SBATCH --partition=YOUR_PARTITION   # set this to your cluster's partition/account
# Define the assemblies directory
ASSEMBLIES_DIR="/path/to/your/analysis_directory/FirstManuscriptTreeComparisons/Kmers_SpeciesGenusFamily/Hypocreaceae/BUSCO_Hypocreaceae/genomes"

# Define the list of same-species accessions (all Trichoderma reesei)
SameSpecies=("GCA_000167675.2" "GCA_002006585.1"
               "GCA_052818895.1" "GCA_047716485.1"
               "GCA_047716445.1" "GCA_047716435.1"
               "GCA_047716415.1" "GCA_042257965.1"
               "GCA_038428105.1" "GCA_028871515.1"
               "GCA_016806875.1" "GCA_028871095.1"
               "GCA_016806815.1" "GCA_028870775.1")

# Create directory
mkdir -p Assemblies_SameSpecies

# Function to create symlinks for a given category
create_symlinks() {
    local category_dir="$1"
    shift
    local accessions=("$@")

    for accession in "${accessions[@]}"; do
        filename="${accession}_genomic.fna"
        source_file="$ASSEMBLIES_DIR/$filename"

        # Check if file exists in the assemblies directory
        if [[ -f "$source_file" ]]; then
            # Create symlink in category directory
            ln -sf "$source_file" "$category_dir/$filename"
            echo "Created symlink: $category_dir/$filename -> $source_file"
        else
            echo "Warning: File $source_file not found"
        fi
    done
}

# Create symlinks
echo "=== Creating symlinks for SameSpecies ==="
create_symlinks "Assemblies_SameSpecies" "${SameSpecies[@]}"

echo "=== Symlink creation complete ==="
echo "Symlinks created in:"
echo "  - Assemblies_SameSpecies/ ($(ls -1 Assemblies_SameSpecies/ 2>/dev/null | wc -l) files)"

# Print summary counts
echo ""
echo "=== Summary ==="
echo "SameSpecies accessions: ${#SameSpecies[@]}"
