from collections import defaultdict
from scipy.stats import spearmanr
from collections import Counter
import numpy as np
import pandas as pd
import sys
import os

# pd.set_option("display.max_columns", None)
# pd.set_option("display.max_rows", 50)
# pd.set_option("display.max_colwidth", 50)
# pd.set_option("display.width", None)

# Define root directory
root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root, "..", "src"))
from default import DATAPATH, CONFIGPATH

# Load pathogen info
pathogen_code = sys.argv[1]
df = pd.read_csv(os.path.join(CONFIGPATH, 'pathogens.csv'))
row = df.loc[df["code"].eq(pathogen_code)]
if row.empty: 
    raise SystemExit(f"Unknown code: {pathogen_code}")
pathogen = row.iloc[0]["pathogen"]

# Define output directory
OUTPUT = os.path.join(root, "..", "output")

print("Step 18: evaluating correlations among datasets")

# Path to correlations
PATH_TO_CORRELATIONS = os.path.join(OUTPUT, "mtuberculosis", "correlations")
STRATEGIES = ['A', 'B', "M"]

# Load probs ref data
PROBS_REF = {}
for strategy in STRATEGIES:
    PROBS_REF[strategy] = {}
    if os.path.exists(os.path.join(PATH_TO_CORRELATIONS, strategy)):
        for dataset in sorted(os.listdir(os.path.join(PATH_TO_CORRELATIONS, strategy))):
            name = dataset.replace("_ref_probs.npz", "")
            probs = np.load(os.path.join(PATH_TO_CORRELATIONS, strategy, dataset))['y_prob_ref']
            PROBS_REF[strategy][name] = probs

# Load ChEMBL data for pathogen
ChEMBL_pathogen = pd.read_csv(os.path.join(OUTPUT, pathogen_code, f"{pathogen_code}_ChEMBL_cleaned_data.csv.gz"), low_memory=False)

# Dict mapping assay_id, activity_type and unit to a set of compound ChEMBL IDs
ASSAY_TO_COMPOUNDS = defaultdict(set)
for assay_id, activity_type, unit, compound_chembl_id in ChEMBL_pathogen[["assay_chembl_id", "activity_type", "unit", "compound_chembl_id"]].values:
    ASSAY_TO_COMPOUNDS[(assay_id, activity_type, unit)].add(compound_chembl_id)
del ChEMBL_pathogen

# Get all compounds from pathogen
compounds = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "compound_counts.csv.gz"))
compounds = set(compounds['compound_chembl_id'])

def calculate_spearman(probs1, probs2):
    return spearmanr(probs1, probs2)

def hit_overlap_chance(probs1, probs2, TOP=100):
    N = len(probs1)
    ind1 = set(np.argsort(probs1)[::-1][:TOP])
    ind2 = set(np.argsort(probs2)[::-1][:TOP])
    m = len(ind1.intersection(ind2))
    expected = TOP * TOP/N  # expected intersection size under chance
    return (m - expected) / (TOP - expected)

# Load data
INDIVIDUAL_SELECTED_LM = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "individual_selected_LM.csv"))
MERGED_SELECTED_LM = pd.read_csv(os.path.join(OUTPUT, pathogen_code, "merged_selected_LM.csv"))
print(f"Total number of datasets: {len(INDIVIDUAL_SELECTED_LM) + len(MERGED_SELECTED_LM)}")

# Filtering for ORGANISM
INDIVIDUAL_SELECTED_LM = INDIVIDUAL_SELECTED_LM[INDIVIDUAL_SELECTED_LM['target_type'] == 'ORGANISM'].reset_index(drop=True)
MERGED_SELECTED_LM = MERGED_SELECTED_LM[MERGED_SELECTED_LM['target_type'] == 'ORGANISM'].reset_index(drop=True)
print(f"Total number of ORGANISM datasets: {len(INDIVIDUAL_SELECTED_LM) + len(MERGED_SELECTED_LM)}")

# Defining final number of compounds in A & B datasets (max among qt and mx)
INDIVIDUAL_SELECTED_LM['cpds'] = INDIVIDUAL_SELECTED_LM[['cpds_qt', 'cpds_mx']].max(axis=1)
INDIVIDUAL_SELECTED_LM['positives'] = INDIVIDUAL_SELECTED_LM[['pos_qt', 'pos_mx']].max(axis=1)

# Define names in A & B datasets
names_AB = []
for assay_id, activity_type, unit, dataset_type, cutoff in INDIVIDUAL_SELECTED_LM[["assay_id", "activity_type", "unit", "dataset_type", "cutoff"]].values:
    dty = "qt" if dataset_type == "quantitative" else "mx"
    names_AB.append(f"{assay_id}_{activity_type}_{unit}_{dty}_{cutoff}")
INDIVIDUAL_SELECTED_LM['name'] = names_AB

# Join tables
INDIVIDUAL_SELECTED_LM_tmp = INDIVIDUAL_SELECTED_LM.copy()
INDIVIDUAL_SELECTED_LM_tmp['assay_keys'] = [str(tuple(i)) for i in INDIVIDUAL_SELECTED_LM_tmp[['assay_id', 'activity_type', 'unit']].values]
INDIVIDUAL_SELECTED_LM_tmp = INDIVIDUAL_SELECTED_LM_tmp.drop(columns=['assay_id', 'is_mid_cutoff', 'pos_qt', 'ratio_qt', 'cpds_qt', 'pos_ql', 'ratio_ql', 'cpds_ql', 'pos_mx', 'ratio_mx', 'cpds_mx', 'overlap_mx'])
INDIVIDUAL_SELECTED_LM_tmp['n_assays'] = 1
MERGED_SELECTED_LM_tmp = MERGED_SELECTED_LM.copy()
MERGED_SELECTED_LM_tmp = MERGED_SELECTED_LM_tmp.drop(columns=['direction', 'assay_type', 'bao_label', 'is_mid_cutoff', 'ratio', 'avg', 'std', 'strain', 'target_chembl_id'])
MERGED_SELECTED_LM_tmp = MERGED_SELECTED_LM_tmp.rename(columns={'n_cpds_union': 'cpds'})
MERGED_SELECTED_LM_tmp['label'] = 'M'
MERGED_SELECTED_LM_tmp['dataset_type'] = np.nan
FINAL_DATASETS = pd.concat([INDIVIDUAL_SELECTED_LM_tmp, MERGED_SELECTED_LM_tmp], ignore_index=True)

# Sorting by compounds
FINAL_DATASETS = FINAL_DATASETS.sort_values(by=["label", "cpds"], ascending=[True, False])

# Mapping names to assay keys
name_to_assaykeys = {i: [eval(k) for k in j.split(";")] for i,j in zip(FINAL_DATASETS['name'], FINAL_DATASETS['assay_keys'])}

# Get names
NAMES = FINAL_DATASETS[['label', 'name']].values.tolist()

# Map name to compounds
NAME_TO_COMPOUNDS = {}
for st, name in NAMES:
    assays = name_to_assaykeys[name]
    cpds = set([cpd for assay in assays for cpd in ASSAY_TO_COMPOUNDS[assay]])
    NAME_TO_COMPOUNDS[name] = cpds

RESULTS = []
RESULTS_DICT = dict()
for n1 in NAMES:
    for n2 in NAMES:
        st1, name1 = n1
        st2, name2 = n2
        probs1 = PROBS_REF[st1][name1]
        probs2 = PROBS_REF[st2][name2]
        cpds1 = NAME_TO_COMPOUNDS[name1]
        cpds2 = NAME_TO_COMPOUNDS[name2]
        a = round(calculate_spearman(probs1, probs2).statistic, 4)
        b = round(hit_overlap_chance(probs1, probs2, TOP=1000), 4)
        c = round(hit_overlap_chance(probs1, probs2, TOP=100), 4)
        d = round(len(cpds1.intersection(cpds2)) / min(len(cpds1), len(cpds2)), 4)
        RESULTS.append([st1, name1, st2, name2, a, b, c, d])
        RESULTS_DICT[(st1, name1, st2, name2)] = (a, b, c, d)

RESULTS = pd.DataFrame(RESULTS, columns=['strategy_1', 'name_1', 'strategy_2', 'name_2', 'spearman', 'hit_overlap_1000', 'hit_overlap_100', 'compound_overlap'])
RESULTS.to_csv(os.path.join(OUTPUT, pathogen_code, 'dataset_correlations.csv'), index=False)

print("All correlations evaluated. Prioritizing datasets...")

SELECTED = []

for label, name in FINAL_DATASETS[['label', 'name']].values:

    select = True
    for previously_selected in SELECTED:
        st1, name1 = label, name
        st2, name2 = previously_selected
        a, b, c, d = RESULTS_DICT[(st1, name1, st2, name2)]
        if a+b+c/3 > 0.5 and d > 0.5:
            select = False
            break
    
    if select == True:
        SELECTED.append([label, name])

FINAL_DATASETS['selected'] = [i in SELECTED for i in FINAL_DATASETS[['label', 'name']].values.tolist()]
FINAL_DATASETS.to_csv(os.path.join(OUTPUT, pathogen_code, 'final_datasets.csv'), index=False)

print(f"Number of selected final datasets per strategy: {dict(Counter(FINAL_DATASETS[FINAL_DATASETS['selected']]['label']))}")

FINAL_COVERAGE = {i: set() for i in 'ABM'}
for label, assay_keys in FINAL_DATASETS[['label', 'assay_keys']].values:
    for assay_key in assay_keys.split(";"):
        assay_key = eval(assay_key)
        FINAL_COVERAGE[label].update(ASSAY_TO_COMPOUNDS[assay_key])

print(f"Final coverage A: {round(100* len(FINAL_COVERAGE['A']) / len(compounds), 1)}%")
print(f"Final coverage B: {round(100* len(FINAL_COVERAGE['B']) / len(compounds), 1)}%")
print(f"Final coverage M: {round(100* len(FINAL_COVERAGE['M']) / len(compounds), 1)}%")
print(f"Final coverage ALL: {round(100* len(FINAL_COVERAGE['A'].union(FINAL_COVERAGE['B']).union(FINAL_COVERAGE['M'])) / len(compounds), 1)}%")