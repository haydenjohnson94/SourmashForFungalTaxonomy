#!/usr/bin/env python
# coding: utf-8

# # Compute PyANI ANIb-based ANI for evolutionary paths of fungi

# In[ ]:


import os
import numpy as np
import pandas as pd
from pathlib import Path
import subprocess
import re
import time
import random
import logging
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Configuration
GENOME_DIR = Path("/path/to/your/analysis_directory/KingdomToSpeciesANI/ForANI_Computations/OutputFiles")
CPUS = 32
WORKERS = 32  # Number of parallel BLAST+ jobs
DATASETS_TOOL = "/path/to/your/datasets_tool/datasets"
DOWNLOAD_THREADS = 4
MAX_RETRIES = 5
MIN_DELAY = 1.2
MAX_DELAY = 3.5
PYANI_OUTDIR = GENOME_DIR / "pyani_output"
PYANI_EXEC = "average_nucleotide_identity.py"  
RESULTS_DIR = GENOME_DIR.parent / "ANI_Results" 
RESULTS_FILE = RESULTS_DIR / "accumulated_ANI_results.csv"

def initialize_results_file():
    """Initialize the results file with headers if it doesn't exist."""
    RESULTS_DIR.mkdir(exist_ok=True)  # Create directory if needed
    if not RESULTS_FILE.exists():
        columns = [
            'Target', 'Same Species', 'Same Genus', 'Same Family',
            'Same Order', 'Same Class', 'Same Phylum', 'Same Kingdom',
            'Other Kingdom'
        ]
        pd.DataFrame(columns=columns).to_csv(RESULTS_FILE, index=False)

def save_results(target_accession, ani_df):
    """Save ANI results to the persistent CSV file in the exact required format."""
    try:
        # Read existing results or create new DataFrame
        if RESULTS_FILE.exists():
            results_df = pd.read_csv(RESULTS_FILE, index_col=0)
        else:
            columns = [
                'Target Accession', 'Same Species', 'Same Genus', 'Same Family',
                'Same Order', 'Same Class', 'Same Phylum', 'Same Kingdom',
                'Other Kingdom'
            ]
            results_df = pd.DataFrame(columns=columns)
        
        # Get the target name from the original accession_paths
        target_name = accession_paths[accession_paths['Target Accession'] == target_accession].index[0]
        
        # Create a new row with all required columns
        new_row = {
            'Target Accession': 1.0,  # First column is always 1.0 (self-comparison)
            'Same Species': np.nan,
            'Same Genus': np.nan,
            'Same Family': np.nan,
            'Same Order': np.nan,
            'Same Class': np.nan,
            'Same Phylum': np.nan,
            'Same Kingdom': np.nan,
            'Other Kingdom': np.nan
        }
        
        # Fill in ANI values for each group
        for group_name in ani_df['group'].unique():
            group_df = ani_df[ani_df['group'] == group_name]
            if not group_df.empty:
                avg_ani = group_df['reciprocal_avg'].mean()
                # Map group names to column names exactly as in your example
                column_map = {
                    'same_species': 'Same Species',
                    'same_genus': 'Same Genus',
                    'same_family': 'Same Family',
                    'same_order': 'Same Order',
                    'same_class': 'Same Class',
                    'same_phylum': 'Same Phylum',
                    'same_kingdom': 'Same Kingdom',
                    'other_kingdom': 'Other Kingdom'
                }
                new_row[column_map[group_name]] = avg_ani
        
        # Add the row to the DataFrame
        results_df.loc[target_name] = new_row
        
        # Save with the exact format you want
        results_df.to_csv(RESULTS_FILE, index=True)
        logging.info(f"Successfully saved results for {target_name} to {RESULTS_FILE}")
        
    except Exception as e:
        logging.error(f"CRITICAL ERROR saving results for {target_accession}: {str(e)}")
        # Print the traceback for debugging
        import traceback
        traceback.print_exc()

                      
# Read accession paths (produced by ../FurtherFilteringEvolutionaryPaths.ipynb)
accession_paths = pd.read_csv('/path/to/your/analysis_directory/KingdomToSpeciesANI/filtered1_evolutionary_paths_accessions.csv')
accession_paths = accession_paths.set_index('Unnamed: 0')

def clean_directory():
    """Clean up .fna files and pyani output directory."""
    logging.info(f"Cleaning up .fna files and pyani output in {GENOME_DIR}")
    for item in GENOME_DIR.glob("*.fna"):
        item.unlink()
    if PYANI_OUTDIR.exists():
        shutil.rmtree(PYANI_OUTDIR)

def download_genome(accession, retries=MAX_RETRIES):
    """Download a single genome with retries."""
    for attempt in range(retries):
        try:
            cmd = f"{DATASETS_TOOL} download genome accession {accession} --filename {GENOME_DIR}/{accession}.zip --exclude-atypical"
            subprocess.run(cmd, shell=True, check=True, capture_output=True)
            logging.info(f"Downloaded {accession}")
            return True
        except subprocess.CalledProcessError as e:
            logging.warning(f"Failed to download {accession}, attempt {attempt + 1}/{retries}: {e}")
            time.sleep(random.uniform(MIN_DELAY, MAX_DELAY))
    logging.error(f"Failed to download {accession} after {retries} attempts")
    return False

def extract_and_organize(accession):
    """Extract and organize downloaded genome files."""
    zip_path = GENOME_DIR / f"{accession}.zip"
    try:
        subprocess.run(f"unzip -o {zip_path} -d {GENOME_DIR}/{accession}", shell=True, check=True)
        for fna in (GENOME_DIR / accession).glob("**/*_genomic.fna"):
            shutil.move(str(fna), GENOME_DIR / f"{accession}_{fna.stem}.fna")
        shutil.rmtree(GENOME_DIR / accession)
        zip_path.unlink()
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to extract/organize {accession}: {e}")

def run_pyani(fna_files, ordered_accessions, accession_paths):
    """Run all required pairwise ANIb comparisons with exact column matching."""
    # Get target accession and its file
    target_accession = ordered_accessions[0]
    try:
        target_fna = next(f for f in fna_files if target_accession in f.name)
    except StopIteration:
        logging.error(f"Target FNA file missing for {target_accession}")
        return None, None  # Return tuple to match expected return

    # Create fresh output directory structure
    try:
        if PYANI_OUTDIR.exists():
            shutil.rmtree(PYANI_OUTDIR)
        PYANI_OUTDIR.mkdir(parents=True)
        (PYANI_OUTDIR / "intermediate").mkdir()
        (PYANI_OUTDIR / "reciprocal_averages").mkdir()
    except OSError as e:
        logging.error(f"Directory setup failed: {e}")
        return None, None

    # Get the target row from accession_paths
    try:
        target_row = accession_paths.loc[accession_paths['Target Accession'] == target_accession].iloc[0]
    except IndexError:
        logging.error(f"Target accession {target_accession} not found in DataFrame")
        return None, None

    # Define comparison groups
    comparison_groups = {
        'same_species': 'Same Species',
        'same_genus': 'Same Genus', 
        'same_family': 'Same Family',
        'same_order': 'Same Order',
        'same_class': 'Same Class',
        'same_phylum': 'Same Phylum',
        'same_kingdom': 'Same Kingdom',
        'other_kingdom': 'Other Kingdom'
    }

    # Prepare all comparisons
    comparisons = []
    for group_name, col_name in comparison_groups.items():
        try:
            accessions = str(target_row[col_name]).split(',')
            for acc in accessions:
                acc = acc.strip()
                if acc and acc.lower() != 'nan':
                    try:
                        compare_fna = next(f for f in fna_files if acc in f.name)
                        comparisons.append((target_accession, target_fna, acc, compare_fna, group_name))
                    except StopIteration:
                        logging.warning(f"Missing FNA for {acc} in {group_name}, skipping")
                        continue
        except KeyError:
            logging.warning(f"Column {col_name} not found in DataFrame, skipping group")
            continue

    if not comparisons:
        logging.error("No valid comparisons found")
        return None, None

    # Process comparisons in parallel
    with ThreadPoolExecutor(max_workers=CPUS) as executor:
        futures = []
        for target_acc, t_fna, comp_acc, c_fna, group in comparisons:
            futures.append(executor.submit(
                process_pairwise_comparison,
                target_acc, t_fna, comp_acc, c_fna, group
            ))
        
        # Collect results
        reciprocal_results = []
        for future in as_completed(futures):
            try:
                result = future.result()
                if result:
                    reciprocal_results.append(result)
            except Exception as e:
                logging.error(f"Comparison failed: {str(e)}")
                continue

    # Save and return final results
    if reciprocal_results:
        reciprocal_df = pd.DataFrame(reciprocal_results)
        reciprocal_df.to_csv(PYANI_OUTDIR / "reciprocal_averages.csv", index=False)
        return reciprocal_df, {f.name.split('_')[0]: f for f in fna_files}
    
    logging.error("No successful ANI comparisons completed")
    return None, None

def process_pairwise_comparison(target_acc, target_fna, comp_acc, comp_fna, group_name):
    """Process a single pairwise comparison with full logging and result saving."""
    pair_dir = PYANI_OUTDIR / "intermediate" / f"{target_acc}_vs_{comp_acc}"
    pair_dir.mkdir(exist_ok=True)
    
    # Create a single comparison directory with both genomes
    comp_dir = pair_dir / f"{target_fna.stem}_vs_{comp_fna.stem}"
    comp_dir.mkdir(exist_ok=True)
    
    try:
        # Link both genome files
        for f in [target_fna, comp_fna]:
            try:
                os.symlink(f, comp_dir / f.name)
            except FileExistsError:
                pass
        
        # Run ANIb just once (it computes both directions automatically)
        output_dir = comp_dir / "pyani_output"
        cmd = [
            PYANI_EXEC,
            "-i", str(comp_dir),
            "-o", str(output_dir),
            "--method", "ANIb",
            "--workers", "4",  # Number of CPUs per comparison
            "--force"
        ]
        
        subprocess.run(cmd, check=True, capture_output=True)
        
        # Parse the ANIb output
        ani_file = output_dir / "ANIb_percentage_identity.tab"
        if not ani_file.exists():
            logging.error(f"ANIb output not found at {ani_file}")
            return None
            
        # Read the ANIb matrix
        ani_matrix = pd.read_csv(ani_file, sep="\t", index_col=0)
        
        # Get the two genome names from the matrix
        genome1, genome2 = ani_matrix.index[0], ani_matrix.index[1]
        
        # Extract forward and reverse ANI values
        forward_ani = ani_matrix.loc[genome1, genome2]
        reverse_ani = ani_matrix.loc[genome2, genome1]
        
        # Calculate reciprocal average
        reciprocal_avg = (forward_ani + reverse_ani) / 2
        
        result = {
            'target': target_acc,
            'comparison': comp_acc,
            'group': group_name,
            'forward_ANI': forward_ani,
            'reverse_ANI': reverse_ani,
            'reciprocal_avg': reciprocal_avg
        }
        
        # Save the complete ANIb matrix (contains both directions)
        ani_matrix.to_csv(pair_dir / "ANIb_matrix.csv")
        
        # Save reciprocal result
        pd.DataFrame([result]).to_csv(
            PYANI_OUTDIR / "reciprocal_averages" / f"{target_acc}_vs_{comp_acc}_reciprocal.csv",
            index=False
        )
        return result
        
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed {comp_dir.name}: {e.stderr.decode()}")
        return None
    except Exception as e:
        logging.error(f"Error processing {target_acc} vs {comp_acc}: {str(e)}")
        return None

        

def process_row(row):
    """Process a single row of accession paths."""
    logging.info(f"Processing row: {row.name}")
    clean_directory()
    
    # Get ordered accessions
    ordered_accessions = [row['Target Accession']] + [row[col] for col in accession_paths.columns[1:] if pd.notna(row[col])]
    
    # Download genomes
    logging.info(f"Downloading {len(ordered_accessions)} assemblies using {DOWNLOAD_THREADS} threads...")
    with ThreadPoolExecutor(max_workers=DOWNLOAD_THREADS) as executor:
        executor.map(download_genome, ordered_accessions)
    
    logging.info("Extracting assemblies...")
    for accession in ordered_accessions:
        extract_and_organize(accession)
    
    logging.info("Organizing files...")
    fna_files = list(GENOME_DIR.glob("*.fna"))
    logging.info(f"Completed: {len(fna_files)}/{len(ordered_accessions)} assemblies ready in {GENOME_DIR}")
    
    # Check if target accession's .fna file exists
    target_accession = ordered_accessions[0]
    target_fna_exists = any(re.match(rf"^{target_accession}", fna.stem) for fna in fna_files)
    if not target_fna_exists:
        logging.error(f"Target accession {target_accession} not found in extracted .fna files, likely excluded as atypical. Skipping row {row.name}")
        return
    
    # Run PyANI ANIb
    logging.info("Computing PyANI ANIb...")
    result = run_pyani(fna_files, ordered_accessions, accession_paths)
    
    if result is None:
        logging.error("PyANI ANIb analysis failed, skipping row")
    else:
        ani_df, _ = result  # Unpack the tuple (ignore accession_to_file)
        try:
            # Log each comparison result
            logging.info(f"ANI results for {target_accession}:")
            for _, result_row in ani_df.iterrows():
                logging.info(
                    f"Comparison: {result_row['target']} vs {result_row['comparison']} | "
                    f"Group: {result_row['group']} | "
                    f"Forward ANI: {result_row['forward_ANI']:.2f} | "
                    f"Reverse ANI: {result_row['reverse_ANI']:.2f} | "
                    f"Reciprocal: {result_row['reciprocal_avg']:.2f}"
                )
            
            # Save detailed results to file
            results_file = PYANI_OUTDIR / f"{target_accession}_full_results.tsv"
            ani_df.to_csv(results_file, sep='\t', index=False)
            logging.info(f"Full results saved to {results_file}")
            
            # Save to persistent CSV
            save_results(target_accession, ani_df)
            
        except KeyError as e:
            logging.error(f"Error processing results: {e}")
            logging.info(f"Available columns: {ani_df.columns.tolist()}")

# Initialize the results file at the start
initialize_results_file()

# Process each row
for _, row in accession_paths.iterrows():
    process_row(row)


# In[ ]:




