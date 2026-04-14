"""
Step 01 — Prepare config files for manual curation.

Generates two files that require human review before the pipeline
can proceed to step 05 and step 08:

1. data/chembl_processed/01_activity_std_units_converted.csv — (activity_type, unit) pairs with
   counts, where activity_type is harmonized and units are converted via
   ucum_manual.csv. A curator fills in the manual_curation_direction column
   (1 = higher value means more active, -1 = lower value means more active)
   and saves the result as:
   config/activity_std_units_manual_curation.csv
   This is required by step 08.

2. data/chembl_processed/01_harmonized_types_map.csv — maps each harmonized activity_type to
   the count and list of raw standard_type strings that collapse into it.
   Human-readable reference; not consumed by downstream pipeline steps.

Run this script after step 00 and before step 05/08.

Input:  data/chembl_processed/00_activity_std_units.csv
        config/ucum_manual.csv
        config/synonyms.csv
Output: data/chembl_processed/01_activity_std_units_converted.csv  (requires manual curation -> config/activity_std_units_manual_curation.csv)
        data/chembl_processed/01_harmonized_types_map.csv
"""

import pandas as pd
import sys
import os

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH
from pathogen_utils import harmonize


def prepare_harmonized_types_map(activity_std_units):
    """Build a mapping from raw standard_type strings to their harmonized canonical form.

    Harmonization strips underscores, spaces, dots, slashes and uppercases the result,
    collapsing variants like 'IC 50', 'ic_50', 'IC/50' -> 'IC50'.

    Produces data/01_harmonized_types_map.csv as a human-readable reference (not consumed by downstream pipeline steps).
    """
    unique_types = activity_std_units["standard_type"].dropna().unique()
    flat = pd.DataFrame({
        "old_type": unique_types,
        "type": [harmonize(t) for t in unique_types]
    })
    mapping = (
        flat.groupby("type")["old_type"]
        .agg(count="count", old_types=lambda x: " ; ".join(x))
        .reset_index()
        .sort_values("count", ascending=False, ignore_index=True)
    )
    outfile = os.path.join(DATAPATH, "chembl_processed", "01_harmonized_types_map.csv")
    os.makedirs(os.path.dirname(outfile), exist_ok=True)
    mapping.to_csv(outfile, index=False)

    print(f"  {len(flat)} unique standard_type values -> {len(mapping)} harmonized types")
    print(f"  {(mapping['count'] > 1).sum()} harmonized types collapsed multiple raw variants")
    print(f"  Saved -> {outfile}")


def prepare_unit_curation_file(activity_std_units):
    """Map (standard_type, standard_units) pairs to harmonized activity_type + unit.

    Applies the same harmonization as prepare_harmonized_types_map to standard_type,
    maps standard_units to final_unit via ucum_manual.csv, collapses synonyms via
    synonyms.csv, then re-aggregates counts by (activity_type, unit).

    Produces data/chembl_processed/01_activity_std_units_converted.csv for manual curation.
    A curator should fill in the manual_curation_direction column (1 = higher value
    means more active, -1 = lower value means more active) and save the result as
    config/activity_std_units_manual_curation.csv before running step 08.
    """
    ucum = pd.read_csv(os.path.join(CONFIGPATH, "ucum_manual.csv"))
    unit_to_final_unit = dict(zip(ucum["units"], ucum["final_unit"]))

    df = activity_std_units.copy()
    df["activity_type"] = df["standard_type"].apply(harmonize)
    df["unit"] = df["standard_units"].map(unit_to_final_unit)

    synonyms = pd.read_csv(os.path.join(CONFIGPATH, "synonyms.csv"))
    for canonical, syns in zip(synonyms["activity_type"], synonyms["synonyms"]):
        for syn in syns.split(";"):
            df.loc[df["activity_type"] == syn.strip(), "activity_type"] = canonical

    out = (
        df.groupby(["activity_type", "unit"], dropna=False)["count"]
        .sum()
        .reset_index()
        .sort_values("count", ascending=False, ignore_index=True)
    )
    total = out["count"].sum()
    out["cumulative_prop"] = (out["count"].cumsum() / total).round(3)

    outfile = os.path.join(DATAPATH, "chembl_processed","01_activity_std_units_converted.csv")
    out.to_csv(outfile, index=False)

    n_mapped = out["unit"].notna().sum()
    n_total = len(out)
    print(f"  {n_mapped}/{n_total} (activity_type, unit) pairs mapped to a final_unit")
    print(f"  Saved -> {outfile}")
    print("  ACTION REQUIRED: fill in manual_curation_direction and save as:")
    print(f"  {os.path.join(CONFIGPATH, 'activity_std_units_manual_curation.csv')}")

if __name__ == "__main__":
    print("Step 01")
    activity_std_units = pd.read_csv(
        os.path.join(DATAPATH, "chembl_processed", "00_activity_std_units.csv"),
        low_memory=False
    )
    print("Preparing harmonized activity types map...")
    prepare_harmonized_types_map(activity_std_units)
    print("Preparing unit curation file...")
    prepare_unit_curation_file(activity_std_units)
    print("Done.")
