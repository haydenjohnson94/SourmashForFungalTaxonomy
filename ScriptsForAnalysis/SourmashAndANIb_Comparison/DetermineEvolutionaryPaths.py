#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import pandas as pd
import json
from ete3 import NCBITaxa
import numpy as np


# In[2]:


## First load and create a dataframe with fungal accessions
# Load the JSON file (shared with DetermineAccessionsDataframe_SpeciesGenus.py, kept in SpeciesGenusAnalysis/)
with open('SpeciesGenusAnalysis/fungi_genome_summary_Apr14_2025.json') as f:
    data1 = json.load(f)
# Extract accession and organism name from each record
records1 = []
for report1 in data1['reports']:
    records1.append({
        'accession': report1['accession'],
        'organism_name': report1['organism']['organism_name']
    })
# Create DataFrame
df1 = pd.DataFrame(records1)




## Then load and create a dataframe with oomyete accessions
# Load the JSON file (shared with DetermineAccessionsDataframe_SpeciesGenus.py, kept in SpeciesGenusAnalysis/)
with open('SpeciesGenusAnalysis/oomycete_genome_summary_Apr14_2025.json') as f:
    data2 = json.load(f)
# Extract accession and organism name from each record
records2 = []
for report2 in data2['reports']:
    records2.append({
        'accession': report2['accession'],
        'organism_name': report2['organism']['organism_name']
    })
# Create DataFrame
df2 = pd.DataFrame(records2)



## Combine the fungal and oomycete NCBI genome assembly summaries
df = pd.concat([df1, df2])




# Filter out uncertain/misclassified taxa, those with genus name in brackets
df = df[~df['organism_name'].str.startswith('[')]


# Filter out accessions that begin with GCF, as they are repeated (they are RefSeq)
df = df[~df['accession'].str.startswith('GCF')]

# Filter out accessions that begin with aff. to make further operations simpler
df = df[~df['organism_name'].str.startswith('aff')]




## This section is to get the species name alone. It takes the second word in the organism name, except for special cases like sp. and aff. named accessions.
# Create a temporary split version of the names
split_names = df['organism_name'].str.split()
# Condition 1: When second word is 'sp.' - take everything after genus
sp_mask = split_names.str[1] == 'sp.'
df.loc[sp_mask, 'Species'] = split_names[sp_mask].str[1:].str.join(' ')
# Condition 2: When second word is 'aff.' - take everything after genus
aff_mask = split_names.str[1] == 'aff.'
df.loc[aff_mask, 'Species'] = split_names[aff_mask].str[1:].str.join(' ')
# Default case: Just take the second word (species epithet)
default_mask = ~(sp_mask | aff_mask)
df.loc[default_mask, 'Species'] = split_names[default_mask].str[1]
# Clean up any missing values (fallback)
df['Species'] = df['Species'].fillna(df['organism_name'].str.split().str[1])




# Extract first two words as species name
df['Genus'] = df['organism_name'].str.split().str[:1].str.join(' ')


# In[4]:


ncbi = NCBITaxa()

# Create empty lists for ALL taxonomic ranks 
kingdom_list = []
phylum_list = []
#subphylum_list = [] 
class_list = []
order_list = []
family_list = []

for organism in df['organism_name']:
    try:
        taxid = ncbi.get_name_translator([organism])
        
        if taxid:
            lineage = ncbi.get_lineage(list(taxid.values())[0][0])
            ranks = ncbi.get_rank(lineage)
            names = ncbi.get_taxid_translator(lineage)

            # Extract all ranks (now including subphylum)
            kingdom = [names[t] for t in lineage if ranks.get(t) == 'kingdom']
            phylum = [names[t] for t in lineage if ranks.get(t) == 'phylum']
            #subphylum = [names[t] for t in lineage if ranks.get(t) == 'subphylum']  
            class_ = [names[t] for t in lineage if ranks.get(t) == 'class']
            order = [names[t] for t in lineage if ranks.get(t) == 'order']
            family = [names[t] for t in lineage if ranks.get(t) == 'family']

            # Append results (first element if exists)
            kingdom_list.append(kingdom[0] if kingdom else "Unknown")
            phylum_list.append(phylum[0] if phylum else "Unknown")
            #subphylum_list.append(subphylum[0] if subphylum else "Unknown")  
            class_list.append(class_[0] if class_ else "Unknown")
            order_list.append(order[0] if order else "Unknown")
            family_list.append(family[0] if family else "Unknown")
        else:
            # Handle missing taxid for all ranks
            kingdom_list.append("Unknown")
            phylum_list.append("Unknown")
            #subphylum_list.append("Unknown") 
            class_list.append("Unknown")
            order_list.append("Unknown")
            family_list.append("Unknown")
    except Exception as e:
        print(f"Error processing {organism}: {str(e)}")
        # Append "Unknown" for all ranks on error
        kingdom_list.append("Unknown")
        phylum_list.append("Unknown")
        #subphylum_list.append("Unknown")  
        class_list.append("Unknown")
        order_list.append("Unknown")
        family_list.append("Unknown")

# Add ALL columns to DataFrame 
df['Family'] = family_list
df['Order'] = order_list
df['Class'] = class_list
#df['Subphylum'] = subphylum_list 
df['Phylum'] = phylum_list
df['Kingdom'] = kingdom_list




# Replace unknown Kingdom for Oomycota with SAR
df.loc[df['Phylum'] == 'Oomycota', 'Kingdom'] = 'Sar'  



# Replace "Unknown" Class for Oomycetes with Peronosporomycetes, as I believe all are of this class
df.loc[(df['Phylum'] == 'Oomycota') & (df['Class'] == 'Unknown'), 'Class'] = 'Peronosporomycetes'


# In[5]:


## Delete rows with species name column as 'x', so species crosses
df = df[df['Species'] != 'x']


## Make a final column that is the full species name without strain info
df['Name'] = df['Genus']+ ' '+df['Species']


# In[6]:


pd.set_option('display.max_rows', 5)  # Show all rows

df


# In[ ]:





# In[ ]:





# In[ ]:





# In[11]:


def create_evolutionary_paths(df, random_state=42):
    # Initialize dataframes
    name_paths = pd.DataFrame()
    accession_paths = pd.DataFrame()
    rank_paths = pd.DataFrame()  

    # Initialize lists to collect data
    name_data = []
    accession_data = []
    rank_data = [] 
    columns = []
    
    # Define the hierarchy levels
    hierarchy = [
        ('Name', 'Same Species'),
        ('Genus', 'Same Genus'),
        ('Family', 'Same Family'), 
        ('Order', 'Same Order'),
        ('Class', 'Same Class'),
        ('Phylum', 'Same Phylum'),
        ('Kingdom', 'Same Kingdom'),
        ('Kingdom', 'Other Kingdom')
    ]
    
    for species in df['Name'].value_counts()[df['Name'].value_counts() > 1].index:
        species_group = df[df['Name'] == species]
        
        for i in range(len(species_group) - 1):
            target = species_group.iloc[i]
            same_species = species_group.iloc[i+1]
            
            # Initialize paths for this pair
            name_path = {'Target Species': target['organism_name'], 'Same Species': same_species['organism_name']}
            accession_path = {'Target Accession': target['accession'], 'Same Species': same_species['accession']}
            rank_path = {'Target Species': target['organism_name'], 'Same Species': same_species['organism_name']}  
            
            for j in range(1, len(hierarchy)):
                rank, level_name = hierarchy[j]
                current_value = target[rank]
                
                if j < len(hierarchy) - 1:
                    if rank == 'Genus':
                        candidates = df[(df[rank] == current_value) & (df['Name'] != species)]
                    else:
                        lower_rank = hierarchy[j-1][0]
                        candidates = df[(df[rank] == current_value) & (df[lower_rank] != target[lower_rank])]
                else:
                    other_kingdom = 'Sar' if current_value == 'Fungi' else 'Fungi'
                    candidates = df[df['Kingdom'] == other_kingdom]
                
                if len(candidates) > 0:
                    selected = candidates.sample(1, random_state=random_state).iloc[0]
                    name_path[level_name] = selected['organism_name']
                    accession_path[level_name] = selected['accession']
                    
                    # Build the rank relationship string
                    if j < len(hierarchy) - 1:
                        # Get the current and lower ranks of the selected organism
                        current_rank_value = selected[rank]
                        lower_rank_value = selected[hierarchy[j-1][0]]
                        
                        # Get the lower rank of the target
                        target_lower_rank = target[hierarchy[j-1][0]]
                        
                        # Format: "CurrentRankValue/LowerRankValue (TargetLowerRank)"
                        rank_path[level_name] = f"{current_rank_value}/{lower_rank_value} ({target_lower_rank})"
                    else:
                        # For Other Kingdom case, just show the kingdom and phylum
                        rank_path[level_name] = f"{selected['Kingdom']}/{selected['Phylum']} ({target['Phylum']})"
                else:
                    name_path[level_name] = np.nan
                    accession_path[level_name] = np.nan
                    rank_path[level_name] = np.nan

            col_name = f"{target['Name']} {i}"
            name_data.append(name_path)
            accession_data.append(accession_path)
            rank_data.append(rank_path)  
            columns.append(col_name)
    
    # Create all dataframes at the end
    name_paths = pd.DataFrame(name_data, index=columns)
    accession_paths = pd.DataFrame(accession_data, index=columns)
    rank_paths = pd.DataFrame(rank_data, index=columns)
    
    return name_paths, accession_paths, rank_paths  

# Generate all paths
name_paths, accession_paths, rank_paths = create_evolutionary_paths(df)


# In[ ]:


# Save the evolutionary paths
name_paths.to_csv('evolutionary_paths_names.csv')

# Save the accession paths  
accession_paths.to_csv('evolutionary_paths_accessions.csv')

# Save to CSV
rank_paths.to_csv('taxonomic_rank_relationships.csv')

