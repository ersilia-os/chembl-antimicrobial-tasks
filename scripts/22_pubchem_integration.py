import pandas as pd
import os

root = os.path.dirname(os.path.abspath(__file__))

# Read supported pathogen codes from config
pathogens_df = pd.read_csv(os.path.join(root, "..", "config", "pathogens.csv"))
pathogen_codes = pathogens_df["code"].tolist()

print(f"Step 22: Computing compounds per assay for {len(pathogen_codes)} pathogens")

for pathogen_code in pathogen_codes:

    raw_path = os.path.join(root, "..", "output", pathogen_code, "07_chembl_raw_data.csv.gz")
    if not os.path.exists(raw_path):
        print(f"  Skipping {pathogen_code}: 07_chembl_raw_data.csv.gz not found")
        continue

    chembl_raw = pd.read_csv(raw_path, low_memory=False)

    assay_to_cpds = {}
    for assay_chembl_id, compound_chembl_id in zip(chembl_raw["assay_chembl_id"], chembl_raw["compound_chembl_id"]):
        if assay_chembl_id not in assay_to_cpds:
            assay_to_cpds[assay_chembl_id] = set()
        assay_to_cpds[assay_chembl_id].add(compound_chembl_id)

    cpd_counts = sorted([[i, len(assay_to_cpds[i])] for i in assay_to_cpds], key=lambda x: x[1], reverse=True)
    df = pd.DataFrame(cpd_counts, columns=["assay_chembl_id", "n_compounds"])

    out_path = os.path.join(root, "..", "output", pathogen_code, f"compounds_per_assay_{pathogen_code}.csv")
    df.to_csv(out_path, index=False)
    print(f"  {pathogen_code}: {len(df)} assays saved")
