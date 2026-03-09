"""
Step 00 — Export ChEMBL tables and generate curation summaries.

Exports 10 raw tables from a local PostgreSQL ChEMBL installation to CSV
(data/chembl_activities/). Then generates frequency-count summary CSVs of
key activity fields.

IMPORTANT: The summary files produced here (activity_comments.csv,
standard_text.csv, activity_std_units.csv) informing the manual curation process.
After running this step, a curator must create/update:
  config/activity_comments_manual_curation.csv
  config/standard_text_manual_curation.csv
  config/activity_std_units_curated_manual_curation.csv
These config files are required by step 05 (clean_activities).
"""

import psycopg
import os
import sys
import pandas as pd


root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATABASE_NAME, CHEMBL_USR, CHEMBL_PWD, DATAPATH

output_dir = os.path.join(DATAPATH, "chembl_activities")

activity_tables = [
    "activities",
    "assays",
    "compound_structures",
    "molecule_dictionary",
    "activity_stds_lookup",
    "target_dictionary",
    "docs",
    "assay_parameters",
    "bioassay_ontology",
    "source",
]

os.makedirs(output_dir, exist_ok=True)

def export_table(conn, table):
    """Export a ChEMBL PostgreSQL table to CSV.
    compound_structures is exported with selected columns only — the binary
    CTAB molfile column is excluded as it is large and not needed downstream.
    """
    outfile = os.path.join(output_dir, f"{table}.csv")
    print(f"Exporting {table} -> {outfile}")

    with conn.cursor() as cur, open(outfile, "w", encoding="utf-8", newline="") as f:
        if table == "compound_structures":
            with cur.copy(f"COPY (SELECT molregno, standard_inchi, standard_inchi_key, canonical_smiles FROM public.{table}) TO STDOUT WITH (FORMAT csv, HEADER true)") as copy:
                for data in copy:
                    f.write(data.tobytes().decode("utf-8"))
        else:
            # Export full table
            with cur.copy(f"COPY public.{table} TO STDOUT WITH (FORMAT csv, HEADER true)") as copy:
                for data in copy:
                    f.write(data.tobytes().decode("utf-8"))

def get_files_from_db():
    """Connect to the local ChEMBL PostgreSQL database and export all tables
    listed in activity_tables to CSV files in the output directory.
    """
    with psycopg.connect(
        dbname=DATABASE_NAME,
        user=CHEMBL_USR,
        password=CHEMBL_PWD,
        host="localhost",
        port=5432
    ) as conn:
        for tbl in activity_tables:
            export_table(conn, tbl)

def curate_activity_files():
    """Read activities.csv and generate frequency-count summaries for key fields.

    NaN values are treated as empty strings (fillna("")) so the empty-string
    row in each output CSV represents "no value recorded". 
    standard_text_value is expected to be populated for <5% of rows
    (qualitative assays only); the remaining rows appear as the empty-string
    bucket. Outputs are the basis for manual curation of the config/ files
    required by step 05.
    """
    df = pd.read_csv(os.path.join(output_dir, "activities.csv"), low_memory=False)

    # Activity comments
    s = df['activity_comment'].astype("string").str.strip().str.lower().fillna("")
    out = (
        s.value_counts(dropna=False)
         .rename_axis('activity_comment')
         .reset_index(name='count')
    )
    out['rank'] = out['count'].rank(method='dense', ascending=False).astype(int)
    out = out.sort_values('count', ascending=False, ignore_index=True)
    total_count = out['count'].sum()
    out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)
    out.to_csv(os.path.join(output_dir, "activity_comments.csv"), index=False)
    print(f"Activity comments [counts] stored in {os.path.join(output_dir, 'activity_comments.csv')}")

    # Activity type - standard units
    s = df[["standard_type", "standard_units"]].astype("string").fillna("")
    out = (
    s.value_counts(subset=["standard_type", "standard_units"], dropna=False)
      .reset_index(name="count")
      .sort_values("count", ascending=False, ignore_index=True)
    )
    total_count = out['count'].sum()
    out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)
    out.to_csv(os.path.join(output_dir, "activity_std_units.csv"), index=False)
    print(f"Activity type - unit pairs [counts] stored in {os.path.join(output_dir, 'activity_std_units.csv')}")

    # Standard units
    s = df['standard_units'].astype("string").str.strip().fillna("")
    out = (
        s.value_counts(dropna=False)
         .rename_axis('standard_units')
         .reset_index(name='count')
    )
    total_count = out['count'].sum()
    out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)
    out = out.sort_values('count', ascending=False, ignore_index=True)
    out.to_csv(os.path.join(output_dir, "standard_units.csv"), index=False)
    print(f"Units [counts] stored in {os.path.join(output_dir, 'standard_units.csv')}")

    # Standard text
    s = df['standard_text_value'].astype("string").str.strip().fillna("")
    out = (
        s.value_counts(dropna=False)
         .rename_axis('standard_text_value')
         .reset_index(name='count')
    )
    total_count = out['count'].sum()
    out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)
    out = out.sort_values('count', ascending=False, ignore_index=True)
    out.to_csv(os.path.join(output_dir, "standard_text.csv"), index=False)
    print(f"Standard text [counts] stored in {os.path.join(output_dir, 'standard_text.csv')}")

def curate_assay_files():
    """Extract (assay_id, chembl_id, description) from assays.csv.

    Produces a lightweight human-readable assay reference used, for example,
    during LLM-based assay parameter curation in step 11.
    """
    df = pd.read_csv(os.path.join(output_dir, "assays.csv"), low_memory=False)
    s = df[['assay_id','chembl_id','description']]
    s = s.sort_values("assay_id", ascending=True, ignore_index=True)
    s.to_csv(os.path.join(output_dir, "assay_descriptions.csv"), index=False)
    print(f"Assay descriptions stored in {os.path.join(output_dir, 'assay_descriptions.csv')}")

if __name__ == "__main__":
    print("Step 00")
    print("Getting files from ChEMBL DB")
    get_files_from_db()
    print("Curating activity files")
    curate_activity_files()
    print("Curating assay description")
    curate_assay_files()

