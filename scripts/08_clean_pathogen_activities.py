from collections import Counter
import pandas as pd
import numpy as np
import sys
import os

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import load_pathogen, load_assay_metadata, build_assays_info

print("Step 08")

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

# 1. Load raw pathogen data
print(f"Loading raw pathogen data for {pathogen_code}...")
chembl_pathogen = pd.read_csv(
    os.path.join(OUTPUT, "07_chembl_raw_data.csv.gz"), low_memory=False
)
assays_raw = pd.read_csv(os.path.join(OUTPUT, "07_assays_raw.csv"))
print(f"  Activities: {len(chembl_pathogen)}")
print(f"  Compounds: {chembl_pathogen['compound_chembl_id'].nunique()}")
print(f"  Assay triplets: {len(assays_raw)} ({assays_raw['assay_id'].nunique()} unique assays)")

# 2. Remove non-standardized compounds
chembl_pathogen = chembl_pathogen[~chembl_pathogen["smiles"].isna()].reset_index(drop=True)
print(f"After removing compounds with no SMILES: {len(chembl_pathogen)} activities, {chembl_pathogen['compound_chembl_id'].nunique()} compounds")

# 3. Remove activities with no numeric value and no text flag
chembl_pathogen = chembl_pathogen[(~chembl_pathogen['value'].isna()) | 
                                (chembl_pathogen['text_flag'] != 0)].reset_index(drop=True)
print(f"After removing empty activities: {len(chembl_pathogen)} activities, {chembl_pathogen['compound_chembl_id'].nunique()} compounds")

# 4. Keep only consensus units (or NaN units)
consensus_units = set(pd.read_csv(os.path.join(DATAPATH, "chembl_processed", "01_activity_std_units_converted.csv"))["unit"].dropna())
chembl_pathogen = chembl_pathogen[
    chembl_pathogen["unit"].isin(consensus_units) | chembl_pathogen["unit"].isna()
].reset_index(drop=True)
print(f"After filtering non-consensus units: {len(chembl_pathogen)} activities, {chembl_pathogen['compound_chembl_id'].nunique()} compounds")

# 5. Assign biological direction per (activity_type, unit)
directions = pd.read_csv(os.path.join(CONFIGPATH, 'activity_std_units_manual_curation.csv'))
directions = {(i,j): k for i,j,k in zip(directions['activity_type'], directions['unit'], directions['manual_curation_direction']) if np.isnan(k) == False}
chembl_pathogen['direction'] = [directions[(i,j)] if (i,j) in directions else np.nan 
                                for i,j in zip(chembl_pathogen['activity_type'], chembl_pathogen['unit'])]
count_directions = Counter(chembl_pathogen['direction'].fillna('NaN'))
print(f"Directions assigned. Summary: {count_directions}")
print(f"Assigned directions [-1, 0, +1]: {round((count_directions[1] + count_directions[-1] + count_directions[0]) / len(chembl_pathogen) * 100, 1)}%")
print(f"Assigned directions [-1, +1]: {round((count_directions[1] + count_directions[-1]) / len(chembl_pathogen) * 100, 1)}%")

# 6. Remove unmodelable activities — keep only those with a direction OR an active/inactive text flag
chembl_pathogen = chembl_pathogen[(chembl_pathogen['direction'].isin([1, -1]) == True) |
                                    (chembl_pathogen['text_flag'].isin([1, -1]))].reset_index(drop=True)
print(f"Keeping only activities with a direction [-1,+1] OR active/inactive text_flag")
print(f"Number of activities (compounds): {len(chembl_pathogen)} ({len(set(chembl_pathogen['compound_chembl_id']))})")

# Saving cleaned pathogen activities
chembl_pathogen.to_csv(os.path.join(OUTPUT,f"08_chembl_cleaned_data.csv.gz"), index=False)

# 7. Activity type - unit summary with text flag counts and direction
print("Preparing activity-unit-comments report...")
flagged = chembl_pathogen["text_flag"].isin([1, -1])
activity_summary = (
    chembl_pathogen[["activity_type", "unit", "text_flag"]]
    .assign(flagged=flagged.to_numpy())
    .groupby(["activity_type", "unit"], dropna=False)
    .agg(count=("flagged", "size"), comments=("flagged", "sum"))
    .reset_index()
    .sort_values("count", ascending=False, ignore_index=True)
)
total = activity_summary["count"].sum()
activity_summary["cumulative_prop"] = (activity_summary["count"].cumsum() / total).round(3)
activity_summary["direction"] = [
    directions.get((at, u), np.nan)
    for at, u in zip(activity_summary["activity_type"], activity_summary["unit"])
]
activity_summary.to_csv(os.path.join(OUTPUT, "08_activity_type_unit_comments.csv"), index=False)

# 8. Build cleaned assay summary table
print("Cleaning individual assay information...")
pathogen_chemical_space = set(pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"])
assay_to_src_id, assay_to_bao_format, src_id_to_src_short_name, bao_id_to_label = load_assay_metadata()

assays_info = build_assays_info(
    chembl_pathogen, pathogen_chemical_space,
    assay_to_src_id, assay_to_bao_format,
    src_id_to_src_short_name, bao_id_to_label,
    directions=directions,
)
assays_info.to_csv(os.path.join(OUTPUT, "08_assays_cleaned.csv"), index=False)
print(f"  Assay triplets written: {len(assays_info)}")
