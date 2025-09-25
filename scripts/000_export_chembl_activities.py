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
    "activity_properties",
    "activity_stds_lookup",
    "activity_supp",
    "activity_supp_map",
    "activity_smid",
    "action_type",
    "assays"
]

os.makedirs(output_dir, exist_ok=True)


def export_table(conn, table):
    outfile = os.path.join(output_dir, f"{table}.csv")
    print(f"Exporting {table} -> {outfile}")

    with conn.cursor() as cur, open(outfile, "w", encoding="utf-8", newline="") as f:
        # psycopg3 COPY streaming
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

    s = df['activity_comment'].astype("string").str.strip().str.lower().fillna("") #To Discuss, make lower?
    out = (
        s.value_counts(dropna=False)
         .rename_axis('activity_comment')
         .reset_index(name='count')
    )
    out['rank'] = out['count'].rank(method='dense', ascending=False).astype(int)
    out = out.sort_values('count', ascending=False, ignore_index=True)
    out.to_csv(os.path.join(output_dir, "activity_comments.csv"), index=False)

    s = df[["standard_type", "standard_units"]].astype("string").fillna("")
    out = (
    s.value_counts(subset=["standard_type", "standard_units"], dropna=False)
      .reset_index(name="count")
      .sort_values("count", ascending=False, ignore_index=True)
    )
    out.to_csv(os.path.join(output_dir, "activity_std_units.csv"), index=False)

    s = df['standard_text_value'].astype("string").str.strip().fillna("")
    out = (
        s.value_counts(dropna=False)
         .rename_axis('activity_comment')
         .reset_index(name='count')
    )
    out = out.sort_values('count', ascending=False, ignore_index=True)
    out.to_csv(os.path.join(output_dir, "standard_text.csv"), index=False)

def curate_assay_files():
    df = pd.read_csv(os.path.join(output_dir, "assays.csv"), low_memory=False)
    s = df[['assay_id','chembl_id','description']]
    s = s.sort_values("assay_id", ascending=True, ignore_index=True)
    s.to_csv(os.path.join(output_dir, "assay_descriptions.csv"), index=False)

if __name__ == "__main__":
    get_files_from_db()
    curate_activity_files()
    curate_assay_files()

