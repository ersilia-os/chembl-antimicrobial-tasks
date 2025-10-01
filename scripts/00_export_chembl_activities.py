import psycopg
import os
import sys
import pandas as pd


root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATABASE_NAME, CHEMBL_USR, CHEMBL_PWD, CONFIGPATH

output_dir = os.path.join(CONFIGPATH, "chembl_activities")

activity_tables = [
    "activities",
    "assays",
    "compound_structures",
    "molecule_dictionary",
    "activity_stds_lookup",
    "target_dictionary"
]

os.makedirs(output_dir, exist_ok=True)

def export_table(conn, table):
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
    df = pd.read_csv(os.path.join(output_dir, "activities.csv"), low_memory=False)

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

    s = df[["standard_type", "standard_units"]].astype("string").fillna("")
    out = (
    s.value_counts(subset=["standard_type", "standard_units"], dropna=False)
      .reset_index(name="count")
      .sort_values("count", ascending=False, ignore_index=True)
    )
    total_count = out['count'].sum()
    out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)
    out.to_csv(os.path.join(output_dir, "activity_std_units.csv"), index=False)

    s = df['standard_text_value'].astype("string").str.strip().fillna("")
    out = (
        s.value_counts(dropna=False)
         .rename_axis('activity_comment')
         .reset_index(name='count')
    )
    # total_count = out['count'].sum()
    # out['cumulative_prop'] = (out['count'].cumsum() / total_count).round(3)
    out = out.sort_values('count', ascending=False, ignore_index=True)
    out.to_csv(os.path.join(output_dir, "standard_text.csv"), index=False)

def curate_assay_files():
    df = pd.read_csv(os.path.join(output_dir, "assays.csv"), low_memory=False)
    s = df[['assay_id','chembl_id','description']]
    s = s.sort_values("assay_id", ascending=True, ignore_index=True)
    s.to_csv(os.path.join(output_dir, "assay_descriptions.csv"), index=False)

if __name__ == "__main__":
    print("Getting files from ChEMBL DB")
    get_files_from_db()
    print("Curating activity files")
    curate_activity_files()
    print("Curating assay description")
    curate_assay_files()

