from collections import defaultdict
import pandas as pd
import numpy as np
import sys
import os

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH, DECOY_RATIO
from pathogen_utils import load_pathogen, load_expert_cutoffs
from dataset_utils import make_dataset_filename
from model_utils import load_ecfp_all, load_all_gz_csvs_from_zip, KFoldTrain, TrainRF, add_decoys, downsample_negatives

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


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def filter_assays_for_merging(assay_df, activity_type, unit, direction, assay_type,
                               target_type_curated_extra, strain=None, bao_label=None,
                               target_chembl_id=None):
    """Filter assay_df to rows matching all provided metadata fields.

    Non-string `unit` and `assay_strain_curated` are treated as NaN (missing).
    `bao_label` and `assay_strain_curated` are optional: pass None to skip filtering on that column.
    If `target_chembl_id` is provided, also filters on `target_chembl_id_curated` (SINGLE PROTEIN only).
    """
    mask = (
        (assay_df["activity_type"] == activity_type) &
        (assay_df["unit"].eq(unit) if isinstance(unit, str) else assay_df["unit"].isna()) &
        (assay_df["direction"] == direction) &
        (assay_df["assay_type"] == assay_type) &
        (assay_df["target_type_curated_extra"] == target_type_curated_extra)
    )
    if bao_label is not None:
        mask &= (assay_df["bao_label"] == bao_label)
    if strain is not None:
        mask &= (assay_df["assay_strain_curated"].eq(strain) if isinstance(strain, str) else assay_df["assay_strain_curated"].isna())
    if target_chembl_id is not None:
        mask &= (assay_df["target_chembl_id_curated"] == target_chembl_id)
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
    out["assay_keys"] = out["assay_keys"].apply(
        lambda ks: ";".join("|".join("" if (not isinstance(x, str) and pd.isna(x)) else str(x) for x in k) for k in ks)
    )

    cols = [c for c in out.columns if c != "assay_keys"] + ["assay_keys"]
    return out[cols].sort_values("n_cpds_union", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------

path_to_correlations = os.path.join(OUTPUT, "13_correlations")
os.makedirs(os.path.join(path_to_correlations, "M"), exist_ok=True)

# Load and merge assay metadata tables
print("Loading assay metadata...")
assays_cleaned = pd.read_csv(os.path.join(OUTPUT, "08_assays_cleaned.csv"))
assays_parameters = pd.read_csv(os.path.join(OUTPUT, "09_assays_parameters_full.csv"))
assay_data_info = pd.read_csv(os.path.join(OUTPUT, "12_assay_data_info.csv"))
individual_selected_lm = pd.read_csv(os.path.join(OUTPUT, "14_individual_selected_LM.csv"))

accepted_assays = set(tuple(row) for row in individual_selected_lm[keys].values)

curated_cols = [
    "target_type_curated", "assay_organism_curated", "target_name_curated",
    "target_chembl_id_curated", "assay_strain_curated", "atcc_id", "mutations",
    "known_drug_resistances", "culture_media",
]
assays_parameters = assays_parameters[keys + curated_cols].drop_duplicates(subset=keys, keep="last").reset_index(drop=True)
assays_parameters[["target_chembl_id_curated", "assay_strain_curated"]] = \
    assays_parameters[["target_chembl_id_curated", "assay_strain_curated"]].replace("", np.nan)
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
ecfps = load_ecfp_all(os.path.join(DATAPATH, "chembl_processed", "06_chembl_ecfps.h5"))

print("Loading reference set...")
reference_set = pd.read_csv(os.path.join(OUTPUT, "13_reference_set.csv.gz"))["reference_compounds"].tolist()
x_ref = np.array([ecfps[cid] for cid in reference_set if cid in ecfps])

pathogen_compounds = set(pd.read_csv(os.path.join(OUTPUT, "07_compound_counts.csv.gz"))["compound_chembl_id"])
decoys_pool = set(i for i in ecfps if i not in pathogen_compounds)

expert_cutoffs = load_expert_cutoffs(CONFIGPATH)

# Load all individual datasets from zip archives into memory
print("Loading individual datasets...")
dfs_qt = load_all_gz_csvs_from_zip(os.path.join(OUTPUT, "12_datasets", "datasets_qt.zip"))
dfs_mx = load_all_gz_csvs_from_zip(os.path.join(OUTPUT, "12_datasets", "datasets_mx.zip"))
print(f"  Quantitative: {len(dfs_qt)} | Mixed: {len(dfs_mx)}")

# ---------------------------------------------------------------------------
# Identify merge candidates
# ---------------------------------------------------------------------------

print("Identifying assays to merge...")
not_accepted = assays_merged[~assays_merged["accepted_in_individual_lm"]]

# Only quantitative/mixed assays can contribute data to merged models
not_accepted = not_accepted[not_accepted["dataset_type"].isin(["quantitative", "mixed"])].copy()

filtered_organism = not_accepted[not_accepted["target_type_curated_extra"] == "ORGANISM"].copy()

# Strain-known ORGANISM: group by strain (no bao_label — ORGANISM target type is sufficient)
keys_organism_strain = ["activity_type", "unit", "direction", "assay_type", "target_type_curated_extra", "assay_strain_curated"]
filtered_organism_known = filtered_organism[filtered_organism["assay_strain_curated"].notna()].copy()
to_merge_organism_strain = to_merge_unique_cpds(filtered_organism_known, keys_organism_strain, assay_to_compounds)
to_merge_organism_strain = to_merge_organism_strain[
    (to_merge_organism_strain["n_cpds_union"] >= 100) & (to_merge_organism_strain["n_assays"] > 1)
].reset_index(drop=True)

# NaN-strain ORGANISM: group only by activity metadata (no strain, no bao_label)
keys_organism_no_strain = ["activity_type", "unit", "direction", "assay_type", "target_type_curated_extra"]
filtered_organism_nan = filtered_organism[filtered_organism["assay_strain_curated"].isna()].copy()
to_merge_organism_no_strain = to_merge_unique_cpds(filtered_organism_nan, keys_organism_no_strain, assay_to_compounds)
to_merge_organism_no_strain["assay_strain_curated"] = np.nan
to_merge_organism_no_strain = to_merge_organism_no_strain[
    (to_merge_organism_no_strain["n_cpds_union"] >= 100) & (to_merge_organism_no_strain["n_assays"] > 1)
].reset_index(drop=True)

to_merge_organism = pd.concat([to_merge_organism_strain, to_merge_organism_no_strain], ignore_index=True)
to_merge_organism["name"] = [f"M_ORG{i}" for i in range(len(to_merge_organism))]

keys_single_protein = ["activity_type", "unit", "direction", "assay_type", "target_type_curated_extra", "bao_label", "assay_strain_curated", "target_chembl_id_curated"]
filtered_single_protein = not_accepted[not_accepted["target_type_curated_extra"] == "SINGLE PROTEIN"].copy()
to_merge_single_protein = to_merge_unique_cpds(filtered_single_protein, keys_single_protein, assay_to_compounds)
col = to_merge_single_protein["target_chembl_id_curated"]
to_merge_single_protein = to_merge_single_protein[
    (to_merge_single_protein["n_cpds_union"] >= 100) &
    (to_merge_single_protein["n_assays"] > 1) &
    (col.notna()) & (col != "")
].reset_index(drop=True)
to_merge_single_protein["name"] = [f"M_SP{i}" for i in range(len(to_merge_single_protein))]

print(f"  ORGANISM groups: {len(to_merge_organism)}")
print(f"  SINGLE PROTEIN groups: {len(to_merge_single_protein)}")

# ---------------------------------------------------------------------------
# Merge and model
# ---------------------------------------------------------------------------
# Track merging analysis for detailed comments in step 18
# ---------------------------------------------------------------------------

def parse_assay_key(assay_key_str):
    """Parse assay key string back to tuple."""
    parts = assay_key_str.split("|")
    # Handle assay_id as string, handle NaN units
    assay_id = parts[0]
    activity_type = parts[1]
    unit = None if parts[2] == "" else parts[2]
    return (assay_id, activity_type, unit)

# Initialize merging analysis tracking
merging_analysis = []

# Track all groups attempted (before filtering) for both target types
def track_all_groups(df_all, group_keys, target_type_name, assay_to_compounds):
    """Track all groups and their failure reasons."""
    groups_all = to_merge_unique_cpds(df_all, group_keys, assay_to_compounds)

    for _, group in groups_all.iterrows():
        assay_keys = [parse_assay_key(ak) for ak in group["assay_keys"].split(";")]
        n_assays = group["n_assays"]
        n_cpds_union = group["n_cpds_union"]

        # Determine failure reason
        if n_assays < 2:
            reason = "insufficient_compatible_assays"
        elif n_cpds_union < 100:
            reason = "insufficient_compounds_after_merging"
        else:
            reason = "group_qualified"  # Will be updated if positives are insufficient

        # Add entry for each assay in this group
        for assay_key in assay_keys:
            merging_analysis.append({
                "assay_id": assay_key[0],
                "activity_type": assay_key[1],
                "unit": assay_key[2],
                "target_type": target_type_name,
                "group_size": n_assays,
                "group_compounds": n_cpds_union,
                "failure_reason": reason,
                "group_keys": str(dict(zip(group_keys, [getattr(group, key) for key in group_keys])))
            })

# Track ORGANISM groups
if len(filtered_organism_known) > 0:
    track_all_groups(filtered_organism_known, keys_organism_strain, "ORGANISM", assay_to_compounds)
if len(filtered_organism_nan) > 0:
    track_all_groups(filtered_organism_nan, keys_organism_no_strain, "ORGANISM", assay_to_compounds)

# Track SINGLE PROTEIN groups
if len(filtered_single_protein) > 0:
    track_all_groups(filtered_single_protein, keys_single_protein, "SINGLE PROTEIN", assay_to_compounds)

# ---------------------------------------------------------------------------

merged_lm = []
label_compounds = {"ORGANISM": set(), "SINGLE PROTEIN": set()}
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
        name = merging.name
        assay_keys = merging.assay_keys
        n_assays = merging.n_assays
        n_cpds_union = merging.n_cpds_union

        if target_type == "ORGANISM":
            # For NaN-strain groups: pass strain=None to skip strain filter
            # For strain-known groups: pass the strain string to filter on it
            raw_strain = merging.assay_strain_curated
            filter_strain = raw_strain if isinstance(raw_strain, str) else None
            target_chembl_id = np.nan
            df = filter_assays_for_merging(
                filtered_assays, activity_type, unit, direction, assay_type,
                target_type_curated_extra, strain=filter_strain,
            )
        else:  # SINGLE PROTEIN
            filter_strain = merging.assay_strain_curated  # string or NaN, both handled in filter function
            target_chembl_id = merging.target_chembl_id_curated
            df = filter_assays_for_merging(
                filtered_assays, activity_type, unit, direction, assay_type,
                target_type_curated_extra, strain=filter_strain,
                bao_label=merging.bao_label,
                target_chembl_id=target_chembl_id,
            )
        df = df[df["dataset_type"].isin(["quantitative", "mixed"])].reset_index(drop=True)

        if len(df) == 0:
            continue

        df_quant = df[df["dataset_type"] == "quantitative"].reset_index(drop=True)
        df_mixed = df[df["dataset_type"] == "mixed"].reset_index(drop=True)

        cutoff_list = expert_cutoffs.get((activity_type, unit, target_type_curated_extra, pathogen_code))
        if not cutoff_list:
            print(f"Warning: Missing expert cutoffs for {activity_type}, {unit}, {target_type_curated_extra}, {pathogen_code}")
            continue
        # ---------------------------------------------------------------------------
        # Two-pass fractional contribution filter (5% of union per pass)
        # Pass 1: qualified assays from the full group
        # Pass 2 (rescue): rejected assays from pass 1, filtered among themselves
        # ---------------------------------------------------------------------------
        MIN_FRACTION = 0.05
        raw_keys = [parse_assay_key(ak) for ak in assay_keys.split(";")]
        union_total = len(set(cid for k in raw_keys for cid in assay_to_compounds.get(k, set())))
        qualified_keys = [k for k in raw_keys if len(assay_to_compounds.get(k, set())) >= MIN_FRACTION * union_total]
        rejected_keys = [k for k in raw_keys if k not in set(qualified_keys)]

        # Rescue pass: apply fractional filter among rejected assays
        if len(rejected_keys) >= 2:
            rescue_union = len(set(cid for k in rejected_keys for cid in assay_to_compounds.get(k, set())))
            rescue_keys = [k for k in rejected_keys if len(assay_to_compounds.get(k, set())) >= MIN_FRACTION * rescue_union]
        else:
            rescue_keys = []
        truly_rejected = [k for k in rejected_keys if k not in set(rescue_keys)]

        # Mark truly rejected assays in merging analysis
        for assay_key in truly_rejected:
            for i, entry in enumerate(merging_analysis):
                if (entry["assay_id"] == assay_key[0] and
                    entry["activity_type"] == assay_key[1] and
                    entry["unit"] == assay_key[2] and
                    entry["failure_reason"] == "group_qualified"):
                    merging_analysis[i]["failure_reason"] = "insufficient_fractional_contribution"

        passes = []
        if len(qualified_keys) >= 2:
            passes.append(("", qualified_keys))
        if len(rescue_keys) >= 2:
            passes.append(("_r", rescue_keys))

        if not passes:
            print(f"  Skipping {name}: no viable pass after fractional filter")
            continue

        for pass_suffix, pass_keys in passes:
            pass_assay_ids = {k[0] for k in pass_keys}
            pass_df_quant = df[
                (df["dataset_type"] == "quantitative") & df["assay_id"].isin(pass_assay_ids)
            ].reset_index(drop=True)
            pass_df_mixed = df[
                (df["dataset_type"] == "mixed") & df["assay_id"].isin(pass_assay_ids)
            ].reset_index(drop=True)
            pass_assay_keys_str = ";".join(
                "|".join("" if (not isinstance(x, str) and pd.isna(x)) else str(x) for x in k)
                for k in pass_keys
            )

            for expert_cutoff in cutoff_list:

                name_ = f"{name}{pass_suffix}_{expert_cutoff}"

                # Load and concatenate quantitative datasets
                data_quant_list = [
                    dfs_qt[make_dataset_filename(a, activity_type, unit, "quantitative", expert_cutoff)].assign(assay_id=a)
                    for a in pass_df_quant["assay_id"]
                    if make_dataset_filename(a, activity_type, unit, "quantitative", expert_cutoff) in dfs_qt
                ]
                data_quant = pd.concat(data_quant_list, ignore_index=True) if data_quant_list else pd.DataFrame()

                # Load and concatenate mixed datasets
                data_mixed_list = [
                    dfs_mx[make_dataset_filename(a, activity_type, unit, "mixed", expert_cutoff)].assign(assay_id=a)
                    for a in pass_df_mixed["assay_id"]
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
                n_real_cpds = len(X)

                ratio_before_decoys = n_positives / len(Y)
                min_cpds = 1000 if ratio_before_decoys < 0.5 else 100

                if n_real_cpds < min_cpds:
                    print(f"  Skipping {name_}: too few compounds ({n_real_cpds} < {min_cpds}, ratio={round(ratio_before_decoys, 2)})")
                    for assay_key in pass_keys:
                        for i, entry in enumerate(merging_analysis):
                            if (entry["assay_id"] == assay_key[0] and
                                entry["activity_type"] == assay_key[1] and
                                entry["unit"] == assay_key[2] and
                                entry["failure_reason"] == "group_qualified"):
                                merging_analysis[i]["failure_reason"] = "insufficient_compounds_after_merging"
                                merging_analysis[i]["n_positives"] = n_positives
                                merging_analysis[i]["group_compounds"] = n_real_cpds
                    continue

                if n_positives <= 50:
                    print(f"  Skipping {name_}: too few positives ({n_positives})")
                    for assay_key in pass_keys:
                        for i, entry in enumerate(merging_analysis):
                            if (entry["assay_id"] == assay_key[0] and
                                entry["activity_type"] == assay_key[1] and
                                entry["unit"] == assay_key[2] and
                                entry["failure_reason"] == "group_qualified"):
                                merging_analysis[i]["failure_reason"] = "insufficient_positives_after_merging"
                                merging_analysis[i]["n_positives"] = n_positives
                                merging_analysis[i]["group_compounds"] = n_real_cpds
                    continue

                print(f"  {name_} | {activity_type} | {unit} | cutoff={expert_cutoff} | strain={filter_strain} | target={target_chembl_id}")
                print(f"    Compounds: {n_real_cpds}, Positives: {n_positives} ({round(100 * n_positives / len(Y), 1)}%)")
                label_compounds[target_type].update(data["compound_chembl_id"].tolist())

                if n_positives / len(Y) > 0.5:
                    X, Y, decoy_ids = add_decoys(X, Y, decoys_pool, ecfps, DECOY_RATIO)
                    print(f"    Added {len(decoy_ids)} decoys")
                    print(f"    After decoys: {len(X)} compounds, {n_positives} positives ({round(100 * n_positives / len(Y), 1)}%)")
                    data = pd.concat([data, pd.DataFrame({"compound_chembl_id": decoy_ids, "bin": 0, "smiles": "decoy"})], ignore_index=True)

                # Downsample negatives for modeling if active ratio < 5%
                active_ratio_model = n_positives / len(Y)
                X, Y = downsample_negatives(X, Y)
                if n_positives / len(Y) != active_ratio_model:
                    print(f"    Downsampled negatives for modeling: ratio {active_ratio_model:.3f} → {n_positives/len(Y):.3f}")

                avg_auroc, std_auroc = KFoldTrain(X, Y)
                print(f"    AUROC: {avg_auroc} ± {std_auroc}")

                # Update merging analysis for successful assays
                for assay_key in pass_keys:
                    for i, entry in enumerate(merging_analysis):
                        if (entry["assay_id"] == assay_key[0] and
                            entry["activity_type"] == assay_key[1] and
                            entry["unit"] == assay_key[2] and
                            entry["failure_reason"] == "group_qualified"):
                            merging_analysis[i]["failure_reason"] = "successfully_merged"
                            merging_analysis[i]["merged_group_name"] = name_
                            merging_analysis[i]["group_compounds"] = n_real_cpds
                            merging_analysis[i]["n_positives"] = n_positives
                            merging_analysis[i]["auroc"] = avg_auroc

                merged_lm.append([
                    name_, activity_type, unit, expert_cutoff, direction, assay_type,
                    target_type_curated_extra, filter_strain, target_chembl_id,
                    len(pass_keys), n_real_cpds, n_positives, round(n_positives / len(Y), 3),
                    avg_auroc, std_auroc, pass_assay_keys_str,
                ])

                # Train final model, save reference predictions and merged dataset
                rf = TrainRF(X, Y)
                y_prob_ref = rf.predict_proba(x_ref)[:, 1]
                np.savez_compressed(os.path.join(path_to_correlations, "M", f"{name_}_ref_probs.npz"), y_prob_ref=y_prob_ref)

                outdir = os.path.join(OUTPUT, "12_datasets", "M")
                os.makedirs(outdir, exist_ok=True)
                data = data[data["smiles"] != "decoy"].reset_index(drop=True)
                data.to_csv(os.path.join(outdir, f"{name_}.csv.gz"), index=False, compression="gzip")

# ---------------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------------

merged_lm_df = pd.DataFrame(merged_lm, columns=[
    "name", "activity_type", "unit", "expert_cutoff", "direction", "assay_type",
    "target_type_curated_extra", "strain", "target_chembl_id",
    "n_assays", "n_cpds_union", "positives", "ratio", "avg", "std", "assay_keys",
])
merged_lm_df.to_csv(os.path.join(OUTPUT, "15_merged_LM.csv"), index=False)

# Save detailed merging analysis for script 18 comments
merging_analysis_df = pd.DataFrame(merging_analysis)
if len(merging_analysis_df) > 0:
    # Fill missing entries for assays that failed to merge
    for col, default in [("n_positives", 0), ("auroc", 0.0), ("merged_group_name", "")]:
        if col not in merging_analysis_df.columns:
            merging_analysis_df[col] = default
        else:
            merging_analysis_df[col] = merging_analysis_df[col].fillna(default)

merging_analysis_df.to_csv(os.path.join(OUTPUT, "15_merging_analysis.csv"), index=False)
print(f"Saved merging analysis for {len(merging_analysis_df)} assay attempts")

for tt, label in [("ORGANISM", "ORG"), ("SINGLE PROTEIN", "SP")]:
    sub = merged_lm_df[merged_lm_df["target_type_curated_extra"] == tt]
    n_assays = len(set(s for row in sub["assay_keys"] for s in row.split(";")))
    cpds = label_compounds[tt]
    print(f"{label} — datasets: {len(sub)}, assays: {n_assays}, coverage: {round(100 * len(cpds) / len(pathogen_compounds), 1)}%")

all_cpds = label_compounds["ORGANISM"] | label_compounds["SINGLE PROTEIN"]
print(f"Overall coverage: {round(100 * len(all_cpds) / len(pathogen_compounds), 1)}%")
