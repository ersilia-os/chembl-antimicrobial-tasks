import pandas as pd
import argparse
import os
import sys

pathogen_code = "calbicans"
data_dir = "../output/calbicans_organism"

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root))
from default_parameters import *

parser = argparse.ArgumentParser(description='Binarize fetched pathogen data')
parser.add_argument('--pathogen_code', type=str, help='Pathogen code')
parser.add_argument('--output_dir', type=str, help='Data directory')

args = parser.parse_args()
data_dir = args.output_dir
pathogen_code = args.pathogen_code

print("Selecting tasks by MOD/DIS")

# Load results
df_MOD = pd.read_csv(os.path.join(data_dir, "014_modelability.csv"))
df_DIS = pd.read_csv(os.path.join(data_dir, "015_distinguishability.csv"))
assert df_DIS['task'].tolist() == df_MOD['task'].tolist()
assert df_DIS['num_pos_samples'].tolist() == df_MOD['num_pos_samples'].tolist()
print(f"Original number of tasks: {len(df_MOD)}")

# Rename columns MOD
new_columns = [f"{col}_MOD" for col in df_MOD.columns]
df_MOD.columns = new_columns

# Rename columns DIS
new_columns = [f"{col}_DIS" for col in df_DIS.columns]
df_DIS.columns = new_columns

# Merge dfs, remove 2 cols and rename the others
merged_df = pd.merge(df_MOD, df_DIS, left_on='task_MOD', right_on='task_DIS', how='inner')
merged_df = merged_df.drop(columns=['task_DIS'])
merged_df = merged_df.rename(columns={'task_MOD': 'task'})
merged_df = merged_df.drop(columns=['num_pos_samples_DIS'])
merged_df = merged_df.rename(columns={'num_pos_samples_MOD': 'num_pos_samples'})

# Non selected tasks will be labeled as 0
merged_df['SELECTED'] = [0] * len(merged_df)

# First filtering condition
merged_df.loc[
    (merged_df['SELECTED'] == 0) & 
    (merged_df['auroc_avg_MOD'] > 0.7) & 
    (merged_df['pos:neg_MOD'] < 1), 
    'SELECTED'
] = 1

# Second filtering condition
merged_df.loc[
    (merged_df['SELECTED'] == 0) & 
    (merged_df['auroc_avg_DIS'] > 0.7) & 
    (merged_df['pos:neg_MOD'] > 1), 
    'SELECTED'
] = 2

# Third filtering condition
merged_df.loc[
    (merged_df['SELECTED'] == 0) & 
    (merged_df['auroc_avg_DIS'] > 0.7) & 
    (merged_df['auroc_avg_MOD'] < 0.7) & 
    (merged_df['num_samples_MOD'] > 1000), 
    'SELECTED'
] = 3

merged_df.to_csv(os.path.join(data_dir, "016_selected_tasks_MOD_DIS.csv"), index=False)
print(f"Number of tasks being selected by MOD/DIS: {len(merged_df[merged_df['SELECTED'] != 0])}")