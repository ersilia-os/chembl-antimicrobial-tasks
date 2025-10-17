from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import sys
import os

# Define root directory
root = "."
sys.path.append(os.path.join(root, "..", "src"))
from default import CONFIGPATH

# Loading ChEMBL preprocessed data
print("Loading ChEMBL preprocessed data...")
ChEMBL = pd.read_csv(os.path.join(root, "..", "config", "chembl_processed", "activities_preprocessed.csv"), low_memory=False)
print(f"Original size: {len(ChEMBL)}")

# Activity type - standard units
s = ChEMBL[["activity_type", "unit"]].astype("string").fillna("")
out = (
s.value_counts(subset=["activity_type", "unit"], dropna=False)
    .reset_index(name="count")
    .sort_values("count", ascending=False, ignore_index=True)
)
total_count = out['count'].sum()
out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)
out['manual_curation'] = ""
out.to_csv(os.path.join(CONFIGPATH, "manual_curation", "activity_std_units_curated_manual_curation_.csv"), index=False)
print(f"Total number of unique activity type - standard unit pairs: {len(out)}")