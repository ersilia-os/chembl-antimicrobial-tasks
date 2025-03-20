import argparse
import pandas as pd
from chemblmltools import chembl_activity_target

# Setup command-line argument parsing
parser = argparse.ArgumentParser(description="Retrieve antimicrobial data from ChEMBL.")

parser.add_argument("--username", required=True, help="Database username")
parser.add_argument("--password", required=True, help="Database password")
parser.add_argument("--db_name", required=True, help="Database name")
parser.add_argument("--organism", default="enterobacter", help="Filter by organism (default: enterobacter)")

# Parse arguments
args = parser.parse_args()

# Fetch data from ChEMBL
try:
    df1 = chembl_activity_target(
        db_user=args.username,
        db_password=args.password,
        db_name=args.db_name,
        organism_contains=args.organism,
    )

    if df1.empty:
        print("No data retrieved. Check the query parameters.")
    else:
        print("\n\nData retrieved successfully! Find the example below: \n\n")
        print(df1.head(5))

except Exception as e:
    print(f"Error retrieving data: {e}")
