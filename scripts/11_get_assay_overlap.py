from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import pickle
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import load_pathogen

# Load pathogen info
pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

print("Step 11")

# Loading pathogen data
print(f"Loading chembl cleaned data for {pathogen_code}...")
chembl_pathogen = pd.read_csv(os.path.join(OUTPUT, "08_chembl_cleaned_data.csv.gz"), low_memory=False)
assays_cleaned = pd.read_csv(os.path.join(OUTPUT, '08_assays_cleaned.csv'))
print(f"Cleaned number of assays: {len(assays_cleaned)}")

# Mapping assay_id - activity_type - unit to a set of compound_chembl_ids
assay_to_compounds = {(assay_id, activity_type, unit, doc_chembl_id): set() for assay_id, activity_type, unit, doc_chembl_id in assays_cleaned[['assay_id', 'activity_type', 'unit', 
                                                                                                                                                "doc_chembl_id"]].values}
print(f"Mapping assay to compounds...")
for assay_id, activity_type, unit, compound_chembl_id, doc_chembl_id in tqdm(chembl_pathogen[['assay_chembl_id', 'activity_type', 'unit', 
                                                                                              'compound_chembl_id', 'doc_chembl_id']].values):
    assay_to_compounds[(assay_id, activity_type, unit, doc_chembl_id)].add(compound_chembl_id)

overlap = []
N = 50
items = [i for i in assay_to_compounds if len(assay_to_compounds[i]) >= 50]
print(f"Number of assays with 50 or more compounds: {len(items)}")

# Restrict to assays with more than N compounds. 
for c, (assay_id_1, activity_type_1, unit_1, doc_1) in tqdm(enumerate(items)):
    cpds_1 = assay_to_compounds[(assay_id_1, activity_type_1, unit_1, doc_1)]
    for assay_id_2, activity_type_2, unit_2, doc_2 in items[c:]:
        cpds_2 = assay_to_compounds[(assay_id_2, activity_type_2, unit_2, doc_2)]
        intersection = len(cpds_1.intersection(cpds_2))
        ratio = round(intersection / min(len(cpds_1), len(cpds_2)), 5)
        overlap.append([assay_id_1, activity_type_1, unit_1, doc_1, assay_id_2, activity_type_2, unit_2, doc_2, len(cpds_1), len(cpds_2), intersection, ratio, doc_1 == doc_2])

# Save results
overlap = pd.DataFrame(overlap, columns=["assay_id_1", "activity_type_1", "unit_1", "doc_1", "assay_id_2", "activity_type_2", "unit_2", "doc_2", 
                                         "cpds_1", "cpds_2", "intersection", "ratio", 'same_doc'])
overlap = overlap.sort_values(by='intersection', ascending=False)
overlap.to_csv(os.path.join(OUTPUT, '11_assays_overlap.csv'), index=False)