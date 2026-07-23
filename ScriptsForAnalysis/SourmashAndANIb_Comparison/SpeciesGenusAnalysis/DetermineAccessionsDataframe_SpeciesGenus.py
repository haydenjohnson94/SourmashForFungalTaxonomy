#!/usr/bin/env python
# coding: utf-8

# In[2]:


import os
import pandas as pd
import json
from ete3 import NCBITaxa
import numpy as np


# In[3]:


## First load and create a dataframe with fungal accessions
# Load the JSON file
with open('fungi_genome_summary_Apr14_2025.json') as f:
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
# Load the JSON file
with open('oomycete_genome_summary_Apr14_2025.json') as f:
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




## This section is to get the species name alone. It takes the second word in the organism name, except for the special case of sp.-named accessions.
## Note: aff.-named organisms are handled separately below by dropping them entirely rather than parsing a species name for them.
# Create a temporary split version of the names
split_names = df['organism_name'].str.split()
# Condition 1: When second word is 'sp.' - take everything after genus
sp_mask = split_names.str[1] == 'sp.'
df.loc[sp_mask, 'Species'] = split_names[sp_mask].str[1:].str.join(' ')
# Default case: Just take the second word (species epithet)
default_mask = ~(sp_mask)
df.loc[default_mask, 'Species'] = split_names[default_mask].str[1]
# Clean up any missing values (fallback)
df['Species'] = df['Species'].fillna(df['organism_name'].str.split().str[1])




# Extract first two words as species name
df['Genus'] = df['organism_name'].str.split().str[:1].str.join(' ')




# 
print(f"Total assemblies: {len(df)}")


# In[6]:


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


# In[7]:


## Delete rows with species name column as 'x', so species crosses
df = df[df['Species'] != 'x']

# Remove rows where 'organism_name' contains 'aff.' OR 'cf.'
df = df[~df['organism_name'].str.contains('aff\.|cf\.', case=False, na=False, regex=True)]

## Make a final column that is the full species name without strain info
df['Name'] = df['Genus']+ ' '+df['Species']


# In[8]:


df


# In[9]:


unique_count = df['Name'].nunique()
print(f"Number of unique names: {unique_count}")


# In[ ]:





# In[17]:


def create_pairwise_comparisons(df):
    # Initialize dataframe with all needed columns
    comparisons_df = pd.DataFrame(columns=[
        'Species1', 'Species2', 
        'Accession1', 'Accession2',
        'Comparison_Type'
    ])
    
    # Get all unique species
    species_counts = df['Name'].value_counts()
    
    # First: All within-species pairs
    for species in species_counts[species_counts > 1].index:
        species_group = df[df['Name'] == species]
        accessions = species_group['accession'].tolist()
        names = species_group['organism_name'].tolist()
        
        # Generate all unique pairwise combinations
        for i in range(len(accessions)):
            for j in range(i+1, len(accessions)):
                comparisons_df = pd.concat([
                    comparisons_df,
                    pd.DataFrame([{
                        'Species1': names[i],
                        'Species2': names[j],
                        'Accession1': accessions[i],
                        'Accession2': accessions[j],
                        'Comparison_Type': 'Within_Species'
                    }])
                ], ignore_index=True)
    
    # Second: All between-species pairs within same genus
    genus_groups = df.groupby('Genus')
    
    for genus, genus_df in genus_groups:
        species_in_genus = genus_df['Name'].unique()
        
        if len(species_in_genus) > 1:
            # Get all accessions and names for each species in this genus
            species_data = {
                species: {
                    'accessions': genus_df[genus_df['Name'] == species]['accession'].tolist(),
                    'names': genus_df[genus_df['Name'] == species]['organism_name'].tolist()
                }
                for species in species_in_genus
            }
            
            # Generate all pairs between different species
            species_list = list(species_data.keys())
            for i in range(len(species_list)):
                for j in range(i+1, len(species_list)):
                    # Create all combinations between these two species
                    for acc1, name1 in zip(species_data[species_list[i]]['accessions'], 
                                         species_data[species_list[i]]['names']):
                        for acc2, name2 in zip(species_data[species_list[j]]['accessions'],
                                            species_data[species_list[j]]['names']):
                            comparisons_df = pd.concat([
                                comparisons_df,
                                pd.DataFrame([{
                                    'Species1': name1,
                                    'Species2': name2,
                                    'Accession1': acc1,
                                    'Accession2': acc2,
                                    'Comparison_Type': 'Within_Genus'
                                }])
                            ], ignore_index=True)
    
    return comparisons_df
    
# Generate dataframes
species_comparisons = create_pairwise_comparisons(df)


# In[14]:


# Save the evolutionary paths
species_comparisons.to_csv('SpeciesGenus_names.csv')


# In[ ]:




