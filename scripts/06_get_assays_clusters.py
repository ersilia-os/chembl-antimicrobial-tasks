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
from default import CONFIGPATH, MIN_ASSAY_SIZE

# Loading ChEMBL preprocessed data
print("Loading ChEMBL preprocessed data...")
ChEMBL = pd.read_csv(os.path.join(root, "..", "config", "chembl_processed", "activities_preprocessed.csv"), low_memory=False)
print(f"Original size: {len(ChEMBL)}")
print("Filtering out nan values...")
ChEMBL = ChEMBL[ChEMBL['value'].isna() == False].reset_index(drop=True)
print(f"Size after filtering nan values: {len(ChEMBL)}")

# List of pathogens to process
pathogens = ["Mycobacterium tuberculosis"]
thrs = [0.4, 0.6, 0.85]

# For each pathogen
for pathogen in pathogens:
    
    # Get assays info
    pathogen_code = str(pathogen.split()[0][0] + pathogen.split()[1]).lower()
    print(f"\n\nFiltering for pathogen: {pathogen_code}...")
    PATH_TO_OUTPUT = os.path.join(root, "..", "output", pathogen_code)
    os.makedirs(PATH_TO_OUTPUT, exist_ok=True)
    ChEMBL_ = ChEMBL[ChEMBL['target_organism'].str.contains(pathogen, case=False, na=False) | 
                    ChEMBL['assay_organism'].str.contains(pathogen, case=False, na=False)].reset_index(drop=True)
    
    print((f"Number of activities for {pathogen}: {len(ChEMBL)}"))
    print(f"Calculating clusters for pathogen: {pathogen_code}...")
    ASSAYS_INFO = pd.read_csv(os.path.join(root, "..", "output", pathogen_code, 'assays.csv'))
    ASSAYS_INFO = ASSAYS_INFO[['assay_id', 'activity_type', 'unit', 'activities', 'cpds']].copy()
    print(f"Original number of assays: {len(ASSAYS_INFO)}")
    CLUSTERS = []

    # For each assay
    for assay_id, activity_type, unit in tqdm(ASSAYS_INFO[['assay_id', 'activity_type', 'unit']].values):

        # If unit is nan
        if pd.isna(unit):
            compounds_assay = ChEMBL_[(ChEMBL_['assay_chembl_id'] == assay_id) & (ChEMBL_['activity_type'] == activity_type) & 
                                     (ChEMBL_['unit'].isna() == True) & (ChEMBL_['canonical_smiles'].isna() == False)]['canonical_smiles'].tolist()
        else:
            compounds_assay = ChEMBL_[(ChEMBL_['assay_chembl_id'] == assay_id) & (ChEMBL_['activity_type'] == activity_type) & 
                                     (ChEMBL_['unit'] == unit) & (ChEMBL_['canonical_smiles'].isna() == False)]['canonical_smiles'].tolist()

        # Get list of unique compounds
        compounds_assay = list(set(compounds_assay))

        # if len(compounds_assay) < 100:
        #     CLUSTERS.append([np.nan for _ in thrs])
        # else:

        # Calculate ECFP4s
        fps = bblean.fps_from_smiles(compounds_assay, kind="ecfp4", pack=True, n_features=2048)

        # Clustering using BitBirch
        thr_to_clusters = {}
        for c, thr in enumerate(thrs):
            bb_tree = bblean.BitBirch(branching_factor=128, threshold=thr, merge_criterion="diameter")
            bb_tree.fit(fps)
            clusters = len(bb_tree.get_cluster_mol_ids())
            thr_to_clusters[thr] = clusters
        CLUSTERS.append([thr_to_clusters[i] for i in thr_to_clusters])

    for c, thr in enumerate(thrs):
        ASSAYS_INFO[f'clusters_{thr}'] = np.array(CLUSTERS)[:, c]

    # Save cluster results
    ASSAYS_INFO.to_csv(os.path.join(root, "..", "output", pathogen_code, 'assays_clusters.csv'), index=False)
