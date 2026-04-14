from collections import Counter
import pandas as pd
import sys
import os

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import load_pathogen, load_manual_assays, load_assay_metadata, build_assays_info

print("Step 07")

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

OUTPATH = os.path.join(root, "..", "output", pathogen_code)
os.makedirs(OUTPATH, exist_ok=True)

# 1. Load ChEMBL preprocessed data
print("Loading ChEMBL preprocessed data...")
chembl = pd.read_csv(
    os.path.join(DATAPATH, "chembl_processed", "05_activities_preprocessed.csv"),
    low_memory=False
)
print(f"  Total activities: {len(chembl)}")

# 2. Filter for pathogen (by organism match or manual assay list)
print(f"Filtering for pathogen: {pathogen}...")
manual_assays = load_manual_assays(pathogen_code)
chembl_pathogen = chembl[
    chembl["target_organism"].str.contains(pathogen, case=False, na=False) |
    chembl["assay_organism"].str.contains(pathogen, case=False, na=False) |
    chembl["assay_chembl_id"].isin(manual_assays)
].reset_index(drop=True)

if chembl_pathogen.empty:
    raise SystemExit(f"No activities found for pathogen '{pathogen}'. Check pathogens.csv.")

print(f"  Activities: {len(chembl_pathogen)}")
print(f"  Unique compounds: {chembl_pathogen['compound_chembl_id'].nunique()}")
chembl_pathogen.to_csv(os.path.join(OUTPATH, f"07_chembl_raw_data.csv.gz"), index=False)

# 3. Target organism frequency table
organism_counts = (
    chembl_pathogen["target_organism"]
    .value_counts()
    .reset_index()
)
organism_counts.to_csv(os.path.join(OUTPATH, "07_target_organism_counts.csv"), index=False)

# 4. Compound metadata and chemical space
compound_info = pd.read_csv(
    os.path.join(DATAPATH, "chembl_processed", "02_compound_info.csv"), low_memory=False
)
ik_dict = dict(zip(compound_info["chembl_id"], compound_info["standard_inchi_key"]))
pair_counts = chembl_pathogen[["compound_chembl_id", "smiles"]].value_counts().reset_index(name="count")
pair_counts["inchikey"] = pair_counts["compound_chembl_id"].map(ik_dict)
pair_counts[["compound_chembl_id", "inchikey", "smiles", "count"]].to_csv(
    os.path.join(OUTPATH, "07_compound_counts.csv.gz"), index=False
)

all_smiles = (
    pair_counts[["smiles"]]
    .drop_duplicates()
    .dropna()
    .reset_index(drop=True)
)
all_smiles.to_csv(os.path.join(OUTPATH, "07_all_smiles.csv"), index=False)

# Define chemical space
pathogen_chemical_space = set(pair_counts["compound_chembl_id"])

# 5. Build assay summary table
print("Collecting individual assay information...")
assay_to_src_id, assay_to_bao_format, src_id_to_src_short_name, bao_id_to_label = load_assay_metadata()

assays_info = build_assays_info(
    chembl_pathogen, pathogen_chemical_space,
    assay_to_src_id, assay_to_bao_format,
    src_id_to_src_short_name, bao_id_to_label,
)
assays_info.to_csv(os.path.join(OUTPATH, "07_assays_raw.csv"), index=False)