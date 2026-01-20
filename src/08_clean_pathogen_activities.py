from collections import defaultdict
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

print("Step 08")

# Define output directory
OUTPUT = os.path.join(root, "..", "output")

def only_one(values, name):
    if len(values) != 1:
        raise ValueError(f"Expected exactly one {name}, found {values}")
    return values[0]

# Loading pathogen data
print(f"Loading ChEMBL preprocessed data for {pathogen_code}...")
ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_raw_data.csv.gz"), low_memory=False)
ASSAYS_RAW = pd.read_csv(os.path.join(OUTPUT, pathogen_code, 'assays_raw.csv'))
print(f"Number of activities for {pathogen_code}: {len(ChEMBL_pathogen)}")
print(f"Number of compounds for {pathogen_code}: {len(set(ChEMBL_pathogen['compound_chembl_id']))}")
print(f"Original number of assays: {len(ASSAYS_RAW)} (unique: {len(set(ASSAYS_RAW['assay_id']))})")

