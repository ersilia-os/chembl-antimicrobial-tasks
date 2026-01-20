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

# Load pathogen info
pathogen_code = sys.argv[1]
df = pd.read_csv(os.path.join(CONFIGPATH, 'pathogens.csv'))
row = df.loc[df["code"].eq(pathogen_code)]
if row.empty: 
    raise SystemExit(f"Unknown code: {pathogen_code}")
pathogen = row.iloc[0]["pathogen"]

# Create output directory
OUTPUT = os.path.join(root, "..", "output")

print("Step 10")

# Loading pathogen data
print(f"Loading ChEMBL cleaned data for {pathogen_code}...")
ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
ASSAYS_CLEANED = pd.read_csv(os.path.join(OUTPUT, pathogen_code, 'assays_cleaned.csv'))
print(f"Cleaned number of assays: {len(ASSAYS_CLEANED)}")

# Mapping assay_id - activity_type - unit to a set of compound_chembl_ids
ASSAY_TO_COMPOUNDS = {(assay_id, activity_type, unit, doc_chembl_id): set() for assay_id, activity_type, unit, doc_chembl_id in ASSAYS_CLEANED[['assay_id', 'activity_type', 'unit', 
                                                                                                                                                "doc_chembl_id"]].values}
print(f"Mapping assay to compounds...")
for assay_id, activity_type, unit, compound_chembl_id, doc_chembl_id in tqdm(ChEMBL_pathogen[['assay_chembl_id', 'activity_type', 'unit', 
                                                                                              'compound_chembl_id', 'doc_chembl_id']].values):
    ASSAY_TO_COMPOUNDS[(assay_id, activity_type, unit, doc_chembl_id)].add(compound_chembl_id)

OVERLAP = []
N = 50
items = [i for i in ASSAY_TO_COMPOUNDS if len(ASSAY_TO_COMPOUNDS[i]) >= 50]
print(f"Number of assays with 50 or more compounds: {len(items)}")

# Restrict to assays with more than N compounds. 
for c, (assay_id_1, activity_type_1, unit_1, doc_1) in tqdm(enumerate(items)):
    cpds_1 = ASSAY_TO_COMPOUNDS[(assay_id_1, activity_type_1, unit_1, doc_1)]
    for assay_id_2, activity_type_2, unit_2, doc_2 in items[c:]:
        cpds_2 = ASSAY_TO_COMPOUNDS[(assay_id_2, activity_type_2, unit_2, doc_2)]
        intersection = len(cpds_1.intersection(cpds_2))
        ratio = round(intersection / min(len(cpds_1), len(cpds_2)), 5)
        OVERLAP.append([assay_id_1, activity_type_1, unit_1, doc_1, assay_id_2, activity_type_2, unit_2, doc_2, len(cpds_1), len(cpds_2), intersection, ratio, doc_1 == doc_2])

# Save results
OVERLAP = pd.DataFrame(OVERLAP, columns=["assay_id_1", "activity_type_1", "unit_1", "doc_1", "assay_id_2", "activity_type_2", "unit_2", "doc_2", 
                                         "cpds_1", "cpds_2", "intersection", "ratio", 'same_doc'])
OVERLAP = OVERLAP.sort_values(by='intersection', ascending=False)
OVERLAP.to_csv(os.path.join(OUTPUT, pathogen_code, 'assays_overlap.csv'), index=False)