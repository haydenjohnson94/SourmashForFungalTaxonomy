#!/usr/bin/env python
# coding: utf-8

# # Compute sourmash inferred ANI for evolutionary paths of fungi

# #### Import dependencies

# In[2]:


import sourmash
import subprocess
import tempfile
import pandas as pd
import numpy as np
import os
from multiprocessing import Pool
from pathlib import Path


# In[ ]:





# In[5]:


### 

import os
import numpy as np
import pandas as pd
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import subprocess
from functools import partial
import re
import time
import random
import logging
import shutil

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Configuration
GENOME_DIR = Path("/path/to/your/analysis_directory/KingdomToSpeciesANI/ForANI_Computations/OutputFiles")
CPUS = 16
SCALED = 1000
K_VALUES = [11, 16, 21, 31, 51]
DATASETS_TOOL = "/path/to/your/datasets_tool/datasets"
DOWNLOAD_THREADS = 4
MAX_RETRIES = 5
MIN_DELAY = 1.2
MAX_DELAY = 3.5
SIGNATURE_TIMEOUT = 1800

# Read accession paths (produced by ../FurtherFilteringEvolutionaryPaths.ipynb)
accession_paths = pd.read_csv('/path/to/your/analysis_directory/KingdomToSpeciesANI/filtered1_evolutionary_paths_accessions.csv')
accession_paths = accession_paths.set_index('Unnamed: 0')


def clean_directory():
    """Clean up .fna, .sig files and directories in GENOME_DIR."""
    logging.info(f"Cleaning up .fna, .sig files and directories in {GENOME_DIR}")
    for item in GENOME_DIR.glob("*.fna"):
        item.unlink()
    for item in GENOME_DIR.glob("*.sig"):
        item.unlink()
    for k_dir in GENOME_DIR.glob("k*"):
        if k_dir.is_dir():
            shutil.rmtree(k_dir)

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

def compute_signature(fna_file, k):
    """Compute sourmash signature for a given k-mer size."""
    sig_path = GENOME_DIR / f"k{k}" / f"{fna_file.stem}_k{k}.sig"
    cmd = f"sourmash sketch dna -p k={k},scaled={SCALED} {fna_file} -o {sig_path}"
    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, timeout=SIGNATURE_TIMEOUT)
        logging.info(f"Completed signature for {fna_file} with k={k}")
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as e:
        logging.error(f"Failed to compute signature for {fna_file} with k={k}: {e}")

def compute_ani(k, ordered_accessions):
    """Compute Jaccard and max-containment ANI for all signatures in a k-mer directory."""
    k_dir = GENOME_DIR / f"k{k}"
    sig_list = k_dir / "sig_list.txt"
    
    # Write signature files in the order of ordered_accessions
    with open(sig_list, "w") as f:
        for acc in ordered_accessions:
            sig_files = list(k_dir.glob(f"{acc}*_k{k}.sig"))
            if sig_files:
                f.write(f"{sig_files[0]}\n")
            else:
                logging.warning(f"No signature file found for accession {acc} with k={k}")
    
    # Compute Jaccard ANI
    jaccard_cmd = f"sourmash compare --from-file {sig_list} -k {k} --ani -o {k_dir}/ani_results_jaccard.npy"
    subprocess.run(jaccard_cmd, shell=True, check=True)
    logging.info(f"Computed Jaccard ANI for k={k}")
    
    # Compute max-containment ANI
    containment_cmd = f"sourmash compare --from-file {sig_list} -k {k} --max-containment --ani -o {k_dir}/ani_results_maxcontainment.npy"
    subprocess.run(containment_cmd, shell=True, check=True)
    logging.info(f"Computed max-containment ANI for k={k}")

def process_row(row):
    """Process a single row of the accession dataframe with ordered ANI results."""
    logging.info(f"Processing row: {row.name}")
    clean_directory()
    
    # Define column order and collect valid (non-NaN) accessions in that order
    columns = ['Target Accession', 'Same Species', 'Same Genus', 'Same Family', 
               'Same Order', 'Same Class', 'Same Phylum', 'Same Kingdom', 'Other Kingdom']
    ordered_accessions = [row[col] for col in columns if pd.notna(row[col])]
    logging.info(f"Downloading {len(ordered_accessions)} assemblies using {DOWNLOAD_THREADS} threads...")
    
    # Download genomes
    with ThreadPoolExecutor(max_workers=DOWNLOAD_THREADS) as executor:
        executor.map(download_genome, ordered_accessions)
    
    logging.info("Extracting assemblies...")
    for accession in ordered_accessions:
        extract_and_organize(accession)
    
    # Check if Target Accession was successfully extracted
    target_accession = row['Target Accession']
    fna_files = list(GENOME_DIR.glob(f"{target_accession}*.fna"))
    if not fna_files:
        logging.error(f"Target accession {target_accession} not found in extracted .fna files, likely excluded as atypical. Skipping row {row.name}")
        return
    
    logging.info("Organizing files...")
    fna_files = list(GENOME_DIR.glob("*.fna"))
    logging.info(f"Completed: {len(fna_files)}/{len(ordered_accessions)} assemblies ready in {GENOME_DIR}")
    
    # Create directories for each k-value
    for k in K_VALUES:
        (GENOME_DIR / f"k{k}").mkdir(exist_ok=True)
        logging.info(f"Processing k={k}")
    
    # Compute signatures
    logging.info("Computing signatures...")
    with ThreadPoolExecutor(max_workers=CPUS) as executor:
        for fna_file in fna_files:
            for k in K_VALUES:
                executor.submit(compute_signature, fna_file, k)
    
    # Compute ANIs
    logging.info("Computing ANIs...")
    for k in K_VALUES:
        compute_ani(k, ordered_accessions)
    
    # Get the ordered list of accessions from the row (including NaN values)
    target_accession = row['Target Accession']
    
    # Create a mapping from accession to filename stem
    accession_to_stem = {}
    for fna in fna_files:
        match = re.match(r"^(GCA_\d+\.\d+)", fna.stem)
        if match:
            acc = match.group(1)
            accession_to_stem[acc] = fna.stem
    
    for k in K_VALUES:
        k_dir = GENOME_DIR / f"k{k}"
        
        # Load the ANI matrices
        jaccard_matrix = np.load(k_dir / "ani_results_jaccard.npy")
        containment_matrix = np.load(k_dir / "ani_results_maxcontainment.npy")
        
        # Load Jaccard labels
        jaccard_labels_path = k_dir / "ani_results_jaccard.npy.labels.txt"
        if not jaccard_labels_path.exists():
            logging.error(f"Jaccard labels file not found: {jaccard_labels_path}")
            continue
        with open(jaccard_labels_path) as f:
            jaccard_labels = [line.strip() for line in f]
        
        # Load max-containment labels
        containment_labels_path = k_dir / "ani_results_maxcontainment.npy.labels.txt"
        if not containment_labels_path.exists():
            logging.error(f"Max-containment labels file not found: {containment_labels_path}")
            continue
        with open(containment_labels_path) as f:
            containment_labels = [line.strip() for line in f]
        
        # Create mapping from accession to matrix index for Jaccard
        jaccard_label_to_index = {}
        for i, label in enumerate(jaccard_labels):
            match = re.match(r".*(GCA_\d+\.\d+)", label)
            if match:
                acc = match.group(1)
                jaccard_label_to_index[acc] = i
            else:
                logging.warning(f"Could not parse accession from Jaccard label: {label}")
        
        # Create mapping from accession to matrix index for max-containment
        containment_label_to_index = {}
        for i, label in enumerate(containment_labels):
            match = re.match(r".*(GCA_\d+\.\d+)", label)
            if match:
                acc = match.group(1)
                containment_label_to_index[acc] = i
            else:
                logging.warning(f"Could not parse accession from max-containment label: {label}")
        
        # Get target index
        if target_accession not in jaccard_label_to_index:
            logging.error(f"Target accession {target_accession} not found in Jaccard ANI matrix labels")
            continue
        if target_accession not in containment_label_to_index:
            logging.error(f"Target accession {target_accession} not found in max-containment ANI matrix labels")
            continue
            
        jaccard_target_idx = jaccard_label_to_index[target_accession]
        containment_target_idx = containment_label_to_index[target_accession]
        
        # Prepare results in the exact order of the dataframe columns
        jaccard_results = {col: np.nan for col in columns}
        containment_results = {col: np.nan for col in columns}
        
        # Target always compares to itself as 1.0
        jaccard_results['Target Accession'] = 1.0
        containment_results['Target Accession'] = 1.0
        
        # Fill in the ANI values in the correct order
        for col in columns[1:]:  # Skip 'Target Accession'
            acc = row[col]
            if pd.notna(acc):
                if acc in jaccard_label_to_index:
                    compare_idx = jaccard_label_to_index[acc]
                    jaccard_results[col] = jaccard_matrix[jaccard_target_idx, compare_idx]
                else:
                    logging.warning(f"Accession {acc} for column {col} not found in Jaccard matrix")
                if acc in containment_label_to_index:
                    compare_idx = containment_label_to_index[acc]
                    containment_results[col] = containment_matrix[containment_target_idx, compare_idx]
                else:
                    logging.warning(f"Accession {acc} for column {col} not found in max-containment matrix")
            else:
                logging.debug(f"Accession for column {col} is NaN")
        
        # Create DataFrames with the same index as the original row
        jaccard_df = pd.DataFrame([jaccard_results], index=[row.name])
        containment_df = pd.DataFrame([containment_results], index=[row.name])
        
        logging.info(f"Ordered Jaccard results for k={k}:\n{jaccard_df}")
        logging.info(f"Ordered Max-Containment results for k={k}:\n{containment_df}")
        
        # Save to CSV files (appending if file exists)
        jaccard_csv = GENOME_DIR / f"ani_results_k{k}_jaccard.csv"
        containment_csv = GENOME_DIR / f"ani_results_k{k}_maxcontainment.csv"
        
        jaccard_df.to_csv(jaccard_csv, mode='a', header=not jaccard_csv.exists())
        containment_df.to_csv(containment_csv, mode='a', header=not containment_csv.exists())
        
        logging.info(f"Appended ordered results to {jaccard_csv} and {containment_csv}")

# Process each row
for _, row in accession_paths.iterrows():
    process_row(row)


# In[ ]:




