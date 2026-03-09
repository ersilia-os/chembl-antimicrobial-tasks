from collections import defaultdict
import pandas as pd
import numpy as np
import random
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import load_pathogen, load_expert_cutoffs
from dataset_utils import make_dataset_filename
from model_utils import load_ecfp_all, load_all_gz_csvs_from_zip, KFoldTrain, TrainRF

pathogen_code = sys.argv[1]
pathogen = load_pathogen(pathogen_code)

print("Step 15")

OUTPUT = os.path.join(root, "..", "output", pathogen_code)

# Shared column lists
keys = ["assay_id", "activity_type", "unit"]
columns_data_info = [
    "target_type_curated_extra", "dataset_type",
    "equal", "higher", "lower",
    "min_", "p1", "p25", "p50", "p75", "p99", "max_",
]

decoy_ratio = 0.1  # target active ratio when adding decoys

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def filter_assays_for_merging(assay_df, activity_type, unit, direction, assay_type,
                               target_type_curated_extra, bao_label, strain,
                               target_chembl_id=None):
    """Filter assay_df to rows matching all provided metadata fields.

    Non-string `unit` and `strain` are treated as NaN (missing).
    If `target_chembl_id` is provided, also filters on that column (SINGLE PROTEIN only).
    """
    mask = (
        (assay_df["activity_type"] == activity_type) &
        (assay_df["unit"].eq(unit) if isinstance(unit, str) else assay_df["unit"].isna()) &
        (assay_df["direction"] == direction) &
        (assay_df["assay_type"] == assay_type) &
        (assay_df["target_type_curated_extra"] == target_type_curated_extra) &
        (assay_df["bao_label"] == bao_label) &
        (assay_df["strain"].eq(strain) if isinstance(strain, str) else assay_df["strain"].isna())
    )
    if target_chembl_id is not None:
        mask &= (assay_df["target_chembl_id"] == target_chembl_id)
    return assay_df[mask].reset_index(drop=True)


def to_merge_unique_cpds(df, group_keys, assay_to_compounds):
    """Group assays by group_keys and compute n_assays, n_cpds_union, and assay_keys.

    assay_keys is stored as a semicolon-separated string of tuple representations,
    e.g. "(assay_id, activity_type, unit);...".
    """
    def collect_assay_keys(block):
        return sorted({tuple(r) for r in block.values})

    def union_size(keys_list):
        u = set()
        for k in keys_list:
            u |= assay_to_compounds.get(k, set())
        return len(u)

    out = (
        df.groupby(group_keys, dropna=False)[["assay_id", "activity_type", "unit"]]
        .apply(collect_assay_keys)
        .reset_index(name="assay_keys")
    )
    out["n_assays"] = out["assay_keys"].apply(len)
    out["n_cpds_union"] = out["assay_keys"].apply(union_size)
    out["assay_keys"] = out["assay_keys"].apply(lambda ks: ";".join("|".join(k) for k in ks))

    cols = [c for c in out.columns if c != "assay_keys"] + ["assay_keys"]
    return out[cols].sort_values("n_cpds_union", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

path_to_correlations = os.path.join(OUTPUT, "correlations")
os.makedirs(os.path.join(path_to_correlations, "M"), exist_ok=True)

# Load and merge assay metadata tables
print("Loading assay metadata...")
assays_cleaned = pd.read_csv(os.path.join(OUTPUT, "08_assays_cleaned.csv"))
assays_parameters = pd.read_csv(os.path.join(OUTPUT, "09_assays_parameters.csv"))
assay_data_info = pd.read_csv(os.path.join(OUTPUT, "12_assay_data_info.csv"))
individual_selected_lm = pd.read_csv(os.path.join(OUTPUT, "14_individual_selected_LM.csv"))

accepted_assays = set(tuple(row) for row in individual_selected_lm[keys].values)

assays_merged = assays_cleaned.merge(assays_parameters, on=keys, how="left", validate="1:1")
assays_merged = assays_merged.merge(assay_data_info[keys + columns_data_info], on=keys, how="left", validate="1:1")
assays_merged["accepted_in_individual_lm"] = [tuple(row) in accepted_assays for row in assays_merged[keys].values]

# Build assay → compound set mapping, then free the full table
print("Mapping assays to compounds...")
chembl = pd.read_csv(os.path.join(OUTPUT, "08_chembl_cleaned_data.csv.gz"), low_memory=False)
assay_to_compounds = defaultdict(set)
for assay_id, activity_type, unit, compound_chembl_id in chembl[["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"]].values:
    assay_to_compounds[(assay_id, activity_type, unit)].add(compound_chembl_id)
del chembl

# Load fingerprints, reference set, and define decoy pool
print("Loading ECFPs...")
ecfps = load_ecfp_all(os.path.join(DATAPATH, "chembl_processed", "ChEMBL_ECFPs.h5"))

print("Loading reference set...")
reference_set = pd.read_csv(os.path.join(OUTPUT, "reference_set.csv.gz"))["reference_compounds"].tolist()
x_ref = np.array([ecfps[cid] for cid in reference_set if cid in ecfps])

pathogen_compounds = set(pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"])
decoys_pool = set(i for i in ecfps if i not in pathogen_compounds)

expert_cutoffs = load_expert_cutoffs(CONFIGPATH)

# Load all individual datasets from zip archives into memory
print("Loading individual datasets...")
dfs_qt = load_all_gz_csvs_from_zip(os.path.join(OUTPUT, "datasets", "datasets_qt.zip"))
dfs_mx = load_all_gz_csvs_from_zip(os.path.join(OUTPUT, "datasets", "datasets_mx.zip"))
print(f"  Quantitative: {len(dfs_qt)} | Mixed: {len(dfs_mx)}")

# ---------------------------------------------------------------------------
# Identify merge candidates
# ---------------------------------------------------------------------------

print("Identifying assays to merge...")
not_accepted = assays_merged[~assays_merged["accepted_in_individual_lm"]]

keys_organism = ["activity_type", "unit", "direction", "assay_type", "target_type_curated_extra", "bao_label", "strain"]
filtered_organism = not_accepted[not_accepted["target_type_curated_extra"] == "ORGANISM"].copy()
to_merge_organism = to_merge_unique_cpds(filtered_organism, keys_organism, assay_to_compounds)
to_merge_organism = to_merge_organism[
    (to_merge_organism["n_cpds_union"] > 1000) & (to_merge_organism["n_assays"] > 1)
].reset_index(drop=True)
to_merge_organism["name"] = [f"M_ORG{i}" for i in range(len(to_merge_organism))]

keys_single_protein = ["activity_type", "unit", "direction", "assay_type", "target_type_curated_extra", "bao_label", "strain", "target_chembl_id"]
filtered_single_protein = not_accepted[not_accepted["target_type_curated_extra"] == "SINGLE PROTEIN"].copy()
to_merge_single_protein = to_merge_unique_cpds(filtered_single_protein, keys_single_protein, assay_to_compounds)
to_merge_single_protein = to_merge_single_protein[
    (to_merge_single_protein["n_cpds_union"] > 1000) &
    (to_merge_single_protein["n_assays"] > 1) &
    (to_merge_single_protein["target_chembl_id"].notna())
].reset_index(drop=True)
to_merge_single_protein["name"] = [f"M_SP{i}" for i in range(len(to_merge_single_protein))]

print(f"  ORGANISM groups: {len(to_merge_organism)}")
print(f"  SINGLE PROTEIN groups: {len(to_merge_single_protein)}")

# ---------------------------------------------------------------------------
# Merge and model
# ---------------------------------------------------------------------------

merged_lm = []
merge_candidates = {
    "ORGANISM": (to_merge_organism, filtered_organism),
    "SINGLE PROTEIN": (to_merge_single_protein, filtered_single_protein),
}

for target_type, (to_merge, filtered_assays) in merge_candidates.items():

    print(f"\n{target_type}...")

    for merging in to_merge.itertuples():

        activity_type = merging.activity_type
        unit = merging.unit
        direction = float(merging.direction)
        assay_type = merging.assay_type
        target_type_curated_extra = merging.target_type_curated_extra
        bao_label = merging.bao_label
        strain = merging.strain
        target_chembl_id = merging.target_chembl_id if target_type == "SINGLE PROTEIN" else np.nan
        name = merging.name
        assay_keys = merging.assay_keys
        n_assays = merging.n_assays
        n_cpds_union = merging.n_cpds_union

        # Filter to matching assays with quantitative or mixed data
        df = filter_assays_for_merging(
            filtered_assays, activity_type, unit, direction, assay_type,
            target_type_curated_extra, bao_label, strain,
            target_chembl_id=target_chembl_id if target_type == "SINGLE PROTEIN" else None,
        )
        df = df[df["dataset_type"].isin(["quantitative", "mixed"])].reset_index(drop=True)

        if len(df) == 0:
            continue

        df_quant = df[df["dataset_type"] == "quantitative"].reset_index(drop=True)
        df_mixed = df[df["dataset_type"] == "mixed"].reset_index(drop=True)

        for expert_cutoff in expert_cutoffs[(activity_type, unit, target_type_curated_extra, pathogen_code)]:

            name_ = f"{name}_{expert_cutoff}"

            # Load and concatenate quantitative datasets
            data_quant_list = [
                dfs_qt[make_dataset_filename(a, activity_type, unit, "quantitative", expert_cutoff)].assign(assay_id=a)
                for a in df_quant["assay_id"]
                if make_dataset_filename(a, activity_type, unit, "quantitative", expert_cutoff) in dfs_qt
            ]
            data_quant = pd.concat(data_quant_list, ignore_index=True) if data_quant_list else pd.DataFrame()

            # Load and concatenate mixed datasets
            data_mixed_list = [
                dfs_mx[make_dataset_filename(a, activity_type, unit, "mixed", expert_cutoff)].assign(assay_id=a)
                for a in df_mixed["assay_id"]
                if make_dataset_filename(a, activity_type, unit, "mixed", expert_cutoff) in dfs_mx
            ]
            if data_mixed_list:
                data_mixed = pd.concat(data_mixed_list, ignore_index=True)
                data_mixed_quant = data_mixed[data_mixed["value"].notna()].reset_index(drop=True)
                data_mixed_qual = data_mixed[data_mixed["value"].isna()].reset_index(drop=True)
            else:
                data_mixed_quant, data_mixed_qual = pd.DataFrame(), pd.DataFrame()

            # Merge quantitative portions
            if len(data_quant) > 0 and len(data_mixed_quant) > 0:
                data = pd.concat([data_quant, data_mixed_quant], ignore_index=True)
            elif len(data_quant) > 0:
                data = data_quant
            elif len(data_mixed_quant) > 0:
                data = data_mixed_quant
            else:
                raise ValueError(f"No quantitative data available for merging {name_}")

            # Deduplicate: keep most active measurement per compound
            ascending = direction == -1
            data = (data.sort_values("value", ascending=ascending)
                        .drop_duplicates("compound_chembl_id", keep="first")
                        .reset_index(drop=True))

            # Append qualitative inactives from mixed datasets
            if len(data_mixed_qual) > 0:
                data = (pd.concat([data, data_mixed_qual], ignore_index=True)
                          .drop_duplicates("compound_chembl_id", keep="first")
                          .reset_index(drop=True))

            X = np.array(data["compound_chembl_id"].map(ecfps).tolist())
            Y = np.array(data["bin"].tolist())
            n_positives = int(Y.sum())

            if n_positives <= 50:
                print(f"  Skipping {name_}: too few positives ({n_positives})")
                continue

            print(f"  {name_} | {activity_type} | {unit} | cutoff={expert_cutoff} | strain={strain} | target={target_chembl_id}")
            print(f"    Compounds: {len(X)}, Positives: {n_positives} ({round(100 * n_positives / len(Y), 1)}%)")

            if n_positives / len(Y) > 0.5:
                n_decoys = int(n_positives / decoy_ratio - (len(Y) - 1))
                print(f"    Adding {n_decoys} decoys")
                rng = random.Random(42)
                decoy_ids = rng.sample(list(decoys_pool), n_decoys)
                X_decoys = np.array([ecfps[i] for i in decoy_ids])
                X = np.vstack([X, X_decoys])
                Y = np.concatenate([Y, np.zeros(len(X_decoys), dtype=Y.dtype)])
                print(f"    After decoys: {len(X)} compounds, {n_positives} positives ({round(100 * n_positives / len(Y), 1)}%)")
                data = pd.concat([data, pd.DataFrame({"compound_chembl_id": decoy_ids, "bin": 0, "smiles": "decoy"})], ignore_index=True)

            avg_auroc, std_auroc = KFoldTrain(X, Y)
            print(f"    AUROC: {avg_auroc} ± {std_auroc}")

            merged_lm.append([
                name_, activity_type, unit, expert_cutoff, direction, assay_type,
                target_type_curated_extra, bao_label, strain, target_chembl_id,
                n_assays, n_cpds_union, n_positives, round(n_positives / len(Y), 3),
                avg_auroc, std_auroc, assay_keys,
            ])

            # Train final model, save reference predictions and merged dataset
            rf = TrainRF(X, Y)
            y_prob_ref = rf.predict_proba(x_ref)[:, 1]
            np.savez_compressed(os.path.join(path_to_correlations, "M", f"{name_}_ref_probs.npz"), y_prob_ref=y_prob_ref)

            outdir = os.path.join(OUTPUT, "datasets", "M")
            os.makedirs(outdir, exist_ok=True)
            data.to_csv(os.path.join(outdir, f"{name_}.csv.gz"), index=False, compression="gzip")

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

merged_lm_df = pd.DataFrame(merged_lm, columns=[
    "name", "activity_type", "unit", "expert_cutoff", "direction", "assay_type",
    "target_type_curated_extra", "bao_label", "strain", "target_chembl_id",
    "n_assays", "n_cpds_union", "positives", "ratio", "avg", "std", "assay_keys",
])
merged_lm_df.to_csv(os.path.join(OUTPUT, "15_merged_LM.csv"), index=False)

for tt, label in [("ORGANISM", "ORG"), ("SINGLE PROTEIN", "SP")]:
    sub = merged_lm_df[merged_lm_df["target_type_curated_extra"] == tt]
    assay_strs = set(s for row in sub["assay_keys"] for s in row.split(";"))
    cpds = set(cpd for s in assay_strs for cpd in assay_to_compounds[tuple(s.split("|"))])
    print(f"{label} — datasets: {len(sub)}, assays: {len(assay_strs)}, coverage: {round(100 * len(cpds) / len(pathogen_compounds), 1)}%")

all_assay_strs = set(s for row in merged_lm_df["assay_keys"] for s in row.split(";"))
all_cpds = set(cpd for s in all_assay_strs for cpd in assay_to_compounds[tuple(s.split("|"))])
print(f"Overall coverage: {round(100 * len(all_cpds) / len(pathogen_compounds), 1)}%")
