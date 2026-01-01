from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))

# List of pathogens to process
pathogens = ["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter",
             "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa",
             "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]
pathogens = ["Acinetobacter baumannii", "Mycobacterium tuberculosis", "Klebsiella pneumoniae"]

def get_pathogen_code(pathogen):
    return str(pathogen.split()[0][0] + pathogen.split()[1]).lower() if len(pathogen.split()) > 1 else pathogen.lower()

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

# For each pathogen
for pathogen in pathogens:

    print("\n\n\n")

    # Loading pathogen data
    pathogen_code = get_pathogen_code(pathogen)
    print(f"Loading ChEMBL cleaned data for {pathogen_code}...")
    ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
    print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
    print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
    ASSAYS_CLEANED = pd.read_csv(os.path.join(OUTPUT, pathogen_code, 'assays_cleaned.csv'))
    print(f"Cleaned number of assays: {len(ASSAYS_CLEANED)}")

    # Mapping assay_id - activity_type - unit to a set of compound_chembl_ids
    ASSAY_TO_COMPOUNDS = {(assay_id, activity_type, unit): set() for assay_id, activity_type, unit in ASSAYS_CLEANED[['assay_id', 'activity_type', 'unit']].values}
    print(f"Mapping assay to compounds...")
    for assay_id, activity_type, unit, compound_chembl_id in tqdm(ChEMBL_pathogen[['assay_chembl_id', 'activity_type', 'unit', 'compound_chembl_id']].values):
        ASSAY_TO_COMPOUNDS[(assay_id, activity_type, unit)].add(compound_chembl_id)

    OVERLAP = []
    N = 50
    items = [i for i in ASSAY_TO_COMPOUNDS if len(ASSAY_TO_COMPOUNDS[i]) >= 50]
    print(f"Number of assays with 50 or more compounds: {len(items)}")

    # Restrict to assays with more than 10? compounds. 
    for c, (assay_id_1, activity_type_1, unit_1) in tqdm(enumerate(items)):
        cpds_1 = ASSAY_TO_COMPOUNDS[(assay_id_1, activity_type_1, unit_1)]
        for assay_id_2, activity_type_2, unit_2 in items[c:]:
            cpds_2 = ASSAY_TO_COMPOUNDS[(assay_id_2, activity_type_2, unit_2)]
            intersection = len(cpds_1.intersection(cpds_2))
            ratio = round(intersection / min(len(cpds_1), len(cpds_2)), 5)
            OVERLAP.append([assay_id_1, activity_type_1, unit_1, assay_id_2, activity_type_2, unit_2, len(cpds_1), len(cpds_2), intersection, ratio])

    # Save results
    OVERLAP = pd.DataFrame(OVERLAP, columns=["assay_id_1", "activity_type_1", "unit_1", "assay_id_2", "activity_type_2", "unit_2", "cpds_1", "cpds_2", "intersection", "ratio"])
    OVERLAP = OVERLAP.sort_values(by='intersection', ascending=False)
    OVERLAP.to_csv(os.path.join(OUTPUT, pathogen_code, 'assay_overlap.csv'), index=False)