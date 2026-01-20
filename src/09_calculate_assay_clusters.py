# Use conda env with bitbirch installed
from collections import Counter
from tqdm import tqdm
import pandas as pd
import numpy as np
import bblean
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

# Thresholds - Tanimoto Coefficient
thrs = [0.3, 0.6, 0.85]

print("Step 09")

print(f"Loading cleaned ChEMBL preprocessed data for {pathogen_code}...")
ChEMBL = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)
print(f"Calculating clusters for pathogen: {pathogen_code}...")
ASSAYS_INFO = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, 'assays_cleaned.csv'))
ASSAYS_INFO = ASSAYS_INFO[['assay_id', 'activity_type', 'unit', 'activities', 'cpds']].copy()
print(f"Original number of assays: {len(ASSAYS_INFO)}")
CLUSTERS = []

# For each assay
for assay_id, activity_type, unit in tqdm(ASSAYS_INFO[['assay_id', 'activity_type', 'unit']].values):

    try:

        # If unit is nan
        if pd.isna(unit):
            compounds_assay = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == activity_type) & 
                                    (ChEMBL['unit'].isna() == True)]['canonical_smiles'].tolist()
        else:
            compounds_assay = ChEMBL[(ChEMBL['assay_chembl_id'] == assay_id) & (ChEMBL['activity_type'] == activity_type) & 
                                    (ChEMBL['unit'] == unit)]['canonical_smiles'].tolist()

        # Get list of unique compounds
        compounds_assay = list(set(compounds_assay))

        # Calculate ECFP4s
        fps = bblean.fps_from_smiles(compounds_assay, kind="ecfp4", pack=True, n_features=2048)

        # Clustering using BitBirch
        thr_to_clusters = {}
        for c, thr in enumerate(thrs):
            try:
                bb_tree = bblean.BitBirch(branching_factor=128, threshold=thr, merge_criterion="diameter")
                bb_tree.fit(fps)
                clusters = len(bb_tree.get_cluster_mol_ids())
                thr_to_clusters[thr] = clusters
            except:
                thr_to_clusters[thr] = np.nan
        CLUSTERS.append([thr_to_clusters[thr] for thr in thrs])

    except:

        CLUSTERS.append([np.nan] * len(thrs))

for c, thr in enumerate(thrs):
    ASSAYS_INFO[f'clusters_{thr}'] = np.array(CLUSTERS)[:, c]

# Save cluster results
ASSAYS_INFO.to_csv(os.path.join(root, "..", "output", pathogen_code, 'assays_clusters.csv'), index=False)
