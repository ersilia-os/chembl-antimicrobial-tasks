from collections import defaultdict
import pandas as pd
import os

pathogen_codes = [
    "abaumannii",
    "calbicans",
    "campylobacter",
    "ecoli",
    "efaecium",
    "enterobacter",
    "hpylori",
    "kpneumoniae",
    "mtuberculosis",
    "ngonorrhoeae",
    "paeruginosa",
    "pfalciparum",
    "saureus",
    "smansoni",
    "spneumoniae",
]

for pathogen_code in pathogen_codes:

    # Get assay IDS
    ChEMBL_pathogen = pd.read_csv(os.path.join('..', 'output', pathogen_code, f"{pathogen_code}_chembl_raw_data.csv.gz"), low_memory=False)
    assay_ids = sorted(set(ChEMBL_pathogen['assay_chembl_id']))
    assay_to_cpds = {i: set() for i in assay_ids}

    # Add compounds to each assay
    for assay_chembl_id, compound_chembl_id in zip(ChEMBL_pathogen['assay_chembl_id'], ChEMBL_pathogen['compound_chembl_id']):
        assay_to_cpds[assay_chembl_id].add(compound_chembl_id)

    # Number of compounds per assay id
    df = sorted([[i, len(assay_to_cpds[i])] for i in assay_to_cpds], key=lambda x: x[1], reverse=True)
    df = pd.DataFrame(df, columns=['assay_chembl_id', 'n_compounds'])
    df.to_csv(os.path.join('..', 'output', pathogen_code, f"compounds_per_assay_{pathogen_code}.csv"), index=False)
    df.to_csv(os.path.join('..', 'output', pathogen_code, f"compounds_per_assay_{pathogen_code}.csv"), index=False)