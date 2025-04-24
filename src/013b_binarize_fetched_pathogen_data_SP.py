import pandas as pd
from tqdm import tqdm
import os
import collections
import numpy as np
import argparse
import sys

root = os.path.dirname(os.path.abspath(__file__))
sys.path.append(os.path.join(root))
from default_parameters import *

parser = argparse.ArgumentParser(description='Binarize fetched pathogen data')
parser.add_argument('--pathogen_code', type=str, help='Pathogen code')
parser.add_argument('--output_dir', type=str, help='Data directory')

###############################
# SINGLE PROTEIN TARGET TYPES #
###############################

args = parser.parse_args()
data_dir = args.output_dir
pathogen_code = args.pathogen_code

# Loading the data
df = pd.read_csv(os.path.join(data_dir, pathogen_code, "012_{0}_cleaned.csv".format(pathogen_code)))
print("Considering only single protein target types")
print("Before: {0}".format(df.shape))
df = df[df["target_type"] == "SINGLE PROTEIN"]
print("After: {0}".format(df.shape))
print("Considering only binding assay types")
print("Before: {0}".format(df.shape))
print(df.value_counts("assay_type"))
df_B = df[df["assay_type"] == "B"]
df_F = df[df["assay_type"] == "F"]
print("After - B: {0}".format(df_B.shape))
print("After - F: {0}".format(df_F.shape))
# df.drop(columns=["assay_type"], inplace=True)

# Labels: binding and functional
labels = {'B': df_B, 'F': df_F}
del df


def pchembl_binarizer(df, prefix):

    """
    Description:
    Binarizes compound activity data based on pChEMBL thresholds and percentiles.
    Generates multiple datasets labeling compounds as active (1) or inactive (0)
    for machine learning tasks.

    Input:
    df (pd.DataFrame): Input DataFrame containing 'pchembl_value', 'pchembl_relation',
                    'inchikey', and 'smiles' columns.
    prefix (str): Prefix to label the resulting datasets.

    Output:
    dict: A dictionary of DataFrames with keys indicating the binarization method
        (e.g., 'prefix_pchembl_value_6.0', 'prefix_pchembl_percentile_10').
        Each DataFrame contains 'inchikey', 'smiles', and binary activity labels.
    """
    
    df = df[df["pchembl_value"].notnull()]
    df = df[df["pchembl_relation"].notnull()]
    data = {}
    for pchembl_cutoff in PCHEMBL_CUTOFFS:
        da = df[df["pchembl_value"] >= pchembl_cutoff]
        da = da[da["pchembl_relation"] != "<"]
        if da.shape[0] < MIN_POSITIVES:
            print("Not enough positives for pchembl cutoff {0}, {1}".format(pchembl_cutoff, da.shape[0]))
            continue
        di = df[df["pchembl_value"] < pchembl_cutoff]
        di = di[di["pchembl_relation"] != ">"]
        actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
        inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
        data["{0}_pchembl_value_{1}".format(prefix, pchembl_cutoff)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                             "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                             "pchembl_value_{0}".format(pchembl_cutoff): [1] * len(actives) + [0] * len(inactives)})
    for percentile in PERCENTILES:
        print("Percentile: {0}".format(percentile))
        df = df[df["pchembl_relation"] != ">"]
        df.loc[df["pchembl_relation"] == "<", "pchembl_value"] = 0
        N = df.shape[0]
        n = int(N * percentile / 100)
        if n < MIN_POSITIVES:
            continue
        df = df.sort_values("pchembl_value", ascending=False)
        da = df.head(n)
        di = df.tail(N - n)
        actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
        inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
        data["{0}_pchembl_percentile_{1}".format(prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                             "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                             "pchembl_percentile_{}".format(percentile): [1] * len(actives) + [0] * len(inactives)})
    print("Collected {0} datasets".format(len(data)))
    return data

def percentage_activity_binarizer(df, prefix):

    """
    Description:
    Binarizes compound bioactivity data based on percentage activity thresholds
    and percentiles. Filters and labels compounds as active (1) or inactive (0)
    for classification tasks.

    Input:
    df (pd.DataFrame): DataFrame containing 'standard_value', 'standard_relation',
                    'standard_units', 'direction_flag', 'inchikey', and 'smiles'.
    prefix (str): Prefix to label the resulting binarized datasets.

    Output:
    dict: A dictionary of DataFrames keyed by binarization method
        (e.g., 'prefix_percentage_activity_50', 'prefix_percentage_activity_percentile_10'),
        each with 'inchikey', 'smiles', and binary activity labels.
    """

    df = df[df["standard_value"].notnull()]
    df = df[df["standard_relation"].notnull()]
    df = df[df["standard_units"] == "%"]
    df = df[df["direction_flag"] == 1]
    data = {}
    for percentage_activity_cutoff in PERCENTAGE_ACTIVITY_CUTOFFS:
        da = df[df["standard_value"] >= percentage_activity_cutoff]
        da = da[da["standard_relation"] != "<"]
        if da.shape[0] < MIN_POSITIVES:
            continue
        di = df[df["standard_value"] < percentage_activity_cutoff]
        di = di[di["standard_relation"] != ">"]
        actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
        inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
        data["{0}_percentage_activity_{1}".format(prefix, percentage_activity_cutoff)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                             "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                             "percentage_activity_{0}".format(percentage_activity_cutoff): [1] * len(actives) + [0] * len(inactives)})
    for percentile in PERCENTILES:
        df = df[df["standard_relation"] != ">"]
        df.loc[df["standard_relation"] == "<", "standard_value"] = 0
        N = df.shape[0]
        n = int(N * percentile / 100)
        if n < MIN_POSITIVES:
            continue
        df = df.sort_values("standard_value", ascending=False)
        da = df.head(n)
        di = df.tail(N - n)
        actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
        inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
        data["{0}_percentage_activity_percentile_{1}".format(prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                             "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                             "percentage_activity_percentile_{0}".format(percentile): [1] * len(actives) + [0] * len(inactives)})
    print("Collected {0} datasets".format(len(data)))
    return data

def active_inactive_binarizer(df, prefix):

    """
    Description:
    Creates a binary classification dataset based on labeled activity flags.
    Compounds marked as active (1) or inactive (-1) are retained and labeled
    accordingly for downstream use.

    Input:
    df (pd.DataFrame): DataFrame containing 'activity_flag', 'inchikey', and 'smiles'.
    prefix (str): Prefix used to label the resulting dataset key.

    Output:
    dict: A dictionary with a single DataFrame containing 'inchikey', 'smiles',
        and a binary 'labeled_active' column (1 = active, 0 = inactive).
    """

    df = df[df["activity_flag"] != 0]
    data = {}
    da = df[df["activity_flag"] == 1]
    di = df[df["activity_flag"] == -1]
    actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
    inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
    data["{0}_labeled_active".format(prefix)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                            "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                            "labeled_active": [1] * len(actives) + [0] * len(inactives)})
    print("Collected {0} datasets".format(len(data)))
    return data

def others_binarizer(df, prefix):
    def split_by_percentiles(df_, direction, inner_prefix):
        data = {}
        for percentile in PERCENTILES:
            if direction == 1:
                df_ = df_[df_["standard_relation"] != ">"]
                df_.loc[df_["standard_relation"] == "<", "standard_value"] = 0
                df_ = df_.sort_values("standard_value", ascending=False)
                N = df_.shape[0]
                n = int(N * percentile / 100)
                if n < MIN_POSITIVES:
                    continue
                da = df_.head(n)
                di = df_.tail(N - n)
                actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
                inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
                data["{0}_{1}_percentile_{2}".format(prefix, inner_prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                            "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                            "{0}_percentile_{1}".format(inner_prefix, percentile): [1] * len(actives) + [0] * len(inactives)})
            elif direction == -1:
                df = df_[df_["standard_relation"] != "<"]
                df.loc[df["standard_relation"] == ">", "standard_value"] = 10**15
                df = df.sort_values("standard_value", ascending=True)
                N = df.shape[0]
                n = int(N * percentile / 100)
                if n < MIN_POSITIVES:
                    continue
                da = df.head(n)
                di = df.tail(N - n)
                actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
                inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
                data["{0}_{1}_percentile_{2}".format(prefix, inner_prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in actives] + [x[0] for x in inactives],
                                            "smiles": [x[1] for x in actives] + [x[1] for x in inactives],
                                            "{0}_percentile_{1}".format(inner_prefix, percentile): [1] * len(actives) + [0] * len(inactives)})
            else:
                continue
            return data

    data = {}
    df = df[df["standard_value"].notnull()]
    df = df[df["standard_relation"].notnull()]
    df = df[df["standard_units"].notnull()]
    df = df[df["pchembl_value"].isnull()]
    df = df[df["pchembl_relation"].isnull()]
    df = df[df["standard_units"] != "%"]
    type_units = []
    for t,u in df[["standard_type", "standard_units"]].values:
        type_units += [(t,u)]
    type_units = list(set(type_units))
    for t, u in type_units:
        df_ = df[(df["standard_type"] == t) & (df["standard_units"] == u)]
        if df_.shape[0] < MIN_SIZE_ANY_TASK:
            continue
        direction = list(set(df_["direction_flag"].tolist()))
        if len(direction) != 1:
            continue
        direction = direction[0]
        direction_confidence = list(set(df_["direction_confidence"].tolist()))
        if len(direction_confidence) != 1:
            continue
        direction_confidence = direction_confidence[0]
        if direction_confidence == 0 or direction == 0:
            for k,v in split_by_percentiles(df_, 1, "{0}_{1}_uncertain".format(t.lower(), u.lower())):
                data[k] = v
            for k,v in split_by_percentiles(df_, -1, "{0}_{1}_uncertain".format(t.lower(), u.lower())):
                data[k] = v
        else:
            for k,v in split_by_percentiles(df_, direction, "{0}_{1}".format(t.lower(), u.lower())):
                data[k] = v
    print("Collected {0} datasets".format(len(data)))
    return data

def append_data(all_datasets, data):
    for k,v in data.items():
        if k not in all_datasets:
            all_datasets[k] = v
        else:
            raise Exception("Dataset {0} already exists".format(k))
    return all_datasets

def create_datasets_by_top_assays(df, all_datasets, priority):
    assay_ids = [x for x in df.value_counts("assay_id").index]
    counts = [x for x in df.value_counts("assay_id")]  # Notice this counts the number of rows, not the number of compounds!
    sel_assay_ids = []
    for aid, count in zip(assay_ids, counts):
        # print(aid, count)
        if count >= MIN_SIZE_ASSAY_TASK:
            sel_assay_ids.append(aid)
        else:
            break
    sel_assay_ids = sel_assay_ids[:MAX_NUM_INDEPENDENT_ASSAYS]
    for aid in sel_assay_ids:
        print("Assay ID: {0}".format(aid))
        dt = df[df["assay_id"] == aid]
        activity_types = [x for x in dt.value_counts("standard_type").index]
        counts = [x for x in dt.value_counts("standard_type")]
        sel_activity_types = []
        for at, count in zip(activity_types, counts):
            if count >= MIN_SIZE_ASSAY_SUBTASK:
                sel_activity_types.append(at)
        sel_activity_types = sel_activity_types[:MAX_NUM_ASSAY_SUBTASKS]
        for activity_type in activity_types:
            prefix = "{0}_assay_{1}_{2}".format(priority, aid, activity_type)
            dtt = dt[dt["standard_type"] == activity_type]
            for has_pchembl in [True, False]:
                if has_pchembl:
                    dttp = dtt[dtt["pchembl_value"].notnull()]
                    if dttp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                        continue
                    data = pchembl_binarizer(dttp, prefix=prefix)
                    all_datasets = append_data(all_datasets, data)
                else:
                    dttp = dtt[dtt["pchembl_value"].isnull()]
                    if dttp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                        continue
                    for has_percentage_activity in [True, False]:
                        if has_percentage_activity:
                            dttpp = dttp[dttp["standard_units"] == "%"]
                            if dttpp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                                continue
                            data = percentage_activity_binarizer(dttpp, prefix=prefix)
                            all_datasets = append_data(all_datasets, data)
                        else:
                            dttpp = dttp[dttp["standard_value"].isnull()]
                            if dttpp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                                continue
                            data = others_binarizer(dttpp, prefix=prefix)
                            all_datasets = append_data(all_datasets, data)
    return all_datasets

def create_datasets_by_active_inactive(df, all_datasets, priority):
    all_targetids = sorted(set(df["target_id"].tolist()))
    for targetid in all_targetids:
        prefix = "{0}_all_{1}".format(priority, targetid)
        df_ = df[df["target_id"] == targetid]
        data = active_inactive_binarizer(df_, prefix=prefix)
        all_datasets = append_data(all_datasets, data)
    return all_datasets

def create_datasets_by_major_types(df, all_datasets, priority):
    counter = collections.defaultdict(int)
    for v in df[["target_id", "standard_type", "standard_units"]].values:
        counter[(v[0], v[1].lower(), v[2])] += 1
    selected_units = sorted(counter.items(), key=lambda x: x[1], reverse=True)[:MAX_NUM_INDEPENDENT_ASSAYS]
    print(selected_units)
    # print(selected_units)
    for r in selected_units:
        r = r[0]
        target_id = r[0]
        standard_type = r[1]
        standard_units = r[2]
        dt = df[(df["target_id"] == target_id) & (df["standard_type"].str.lower() == standard_type) & (df["standard_units"] == standard_units)]
        prefix = "{0}_target_{1}_{2}_{3}".format(priority, target_id, standard_type.lower(), standard_units.lower())
        for has_pchembl in [True, False]:
            if has_pchembl:
                dtp = dt[dt["pchembl_value"].notnull()]
                if dtp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                    continue
                data = pchembl_binarizer(dtp, prefix=prefix)
                all_datasets = append_data(all_datasets, data)
            else:
                dtp = dt[dt["pchembl_value"].isnull()]
                if dtp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                    continue
                for has_percentage_activity in [True, False]:
                    if has_percentage_activity:
                        dtpp = dtp[dtp["standard_units"] == "%"]
                        if dtpp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                            continue
                        data = percentage_activity_binarizer(dtpp, prefix=prefix)
                        all_datasets = append_data(all_datasets, data)
                    else:
                        dtpp = dtp[dtp["standard_value"].isnull()]
                        if dtpp.shape[0] < MIN_SIZE_ASSAY_SUBTASK:
                            continue
                        data = others_binarizer(dtpp, prefix=prefix)
                        all_datasets = append_data(all_datasets, data)
    return all_datasets 

def create_datasets_by_all_pchembl_target(df, all_datasets, priority):
    print("Considering only pchembl values")
    all_targetids = sorted(set(df["target_id"].tolist()))
    for targetid in all_targetids:
        print("Target ID pChEMBL: {0}".format(targetid))
        dtp = df[(df["pchembl_value"].notnull()) & (df["target_id"] == targetid)]
        if dtp.shape[0] > MIN_SIZE_ANY_TASK:
            prefix = "{0}_all_{1}".format(priority, targetid)
            data = pchembl_binarizer(dtp, prefix=prefix)
            all_datasets = append_data(all_datasets, data)
    return all_datasets

def create_datasets_by_all_percentage_target(df, all_datasets, priority):
    print("Considering only percentage activity")
    dtp = df[df["standard_units"] == "%"]
    dtp = dtp[dtp["standard_value"].notnull()]
    dtp = dtp[dtp["standard_relation"].notnull()]
    dtp = dtp[dtp["direction_flag"] == 1]
    all_targetids = sorted(set(dtp["target_id"].tolist()))
    for targetid in all_targetids:
        print("Target ID %: {0}".format(targetid))
        dtp_ = dtp[dtp["target_id"] == targetid]
        if dtp_.shape[0] > MIN_SIZE_ANY_TASK:
            prefix = "{0}_all_{1}".format(priority, targetid)
            data = percentage_activity_binarizer(dtp_, prefix=prefix)
            all_datasets = append_data(all_datasets, data)
    return all_datasets

def create_datasets_by_grouping_percentiles_target(df, all_datasets, priority):
    print("Grouping percentiles")
    all_targetids = sorted(set(df["target_id"].tolist()))
    for targetid in all_targetids:
        all_units = []
        data_actives = collections.defaultdict(list)
        data_inactives = collections.defaultdict(list)
        for v in df[["target_id", "standard_type", "standard_units"]].values:
            if v[0] == targetid:
                all_units.append((v[0], v[1].lower(), v[2]))
        all_units = list(set(all_units))
        for r in tqdm(all_units):
            dp = df[(df["target_id"] == r[0]) & (df["standard_type"].str.lower() == r[1]) & (df["standard_units"] == r[2])]
            for direction in [1, -1]:
                if direction == 1:
                    dp = dp[dp["standard_relation"] != ">"]
                    dp.loc[dp["standard_relation"] == "<", "standard_value"] = 0
                    dp = dp.sort_values("standard_value", ascending=False)
                elif direction == -1:
                    dp = dp[dp["standard_relation"] != "<"]
                    dp.loc[dp["standard_relation"] == ">", "standard_value"] = 10**15
                    dp = dp.sort_values("standard_value", ascending=True)
                else:
                    continue
                N = dp.shape[0]
                for percentile in PERCENTILES:
                    n = int(N * percentile / 100)
                    if n == 0: 
                        continue
                    da = dp.head(n)
                    di = dp.tail(N - n)
                    actives = [(ik, smi) for ik, smi in da[["inchikey", "smiles"]].values]
                    inactives = [(ik, smi) for ik, smi in di[["inchikey", "smiles"]].values]
                    data_actives[percentile] += actives
                    data_inactives[percentile] += inactives
        prefix = "{0}_{1}_grouped_percentiles".format(priority, targetid)
        data = {}
        for percentile in PERCENTILES:
            data["{0}_{1}".format(prefix, percentile)] = pd.DataFrame({"inchikey": [x[0] for x in data_actives[percentile]] + [x[0] for x in data_inactives[percentile]],
                                            "smiles": [x[1] for x in data_actives[percentile]] + [x[1] for x in data_inactives[percentile]],
                                            "percentile_{0}".format(percentile): [1] * len(data_actives[percentile]) + [0] * len(data_inactives[percentile])})
        all_datasets = append_data(all_datasets, data)
    return all_datasets


for label in sorted(labels):

    # Binding or Functional
    df = labels[label]

    # Specify tasks directory
    tasks_dir = os.path.join(data_dir, pathogen_code, "013b_raw_tasks_SP", label)
    if not os.path.exists(tasks_dir):
        os.makedirs(tasks_dir)

    all_datasets = {}
    all_datasets = create_datasets_by_top_assays(df, all_datasets, priority=1)
    all_datasets = create_datasets_by_major_types(df, all_datasets, priority=2)
    all_datasets = create_datasets_by_all_pchembl_target(df, all_datasets, priority=3)  # Adapted function, all pchembl values PER TARGET
    all_datasets = create_datasets_by_all_percentage_target(df, all_datasets, priority=4)  # Adapted function, all percentages PER TARGET
    all_datasets = create_datasets_by_grouping_percentiles_target(df, all_datasets, priority=5)  # Adapted function, all percentiles PER TARGET
    all_datasets = create_datasets_by_active_inactive(df, all_datasets, priority=6)

    def disambiguate_data(df):
        ik2smi = {}
        ik2act = collections.defaultdict(list)
        columns = list(df.columns)
        assert len(columns) == 3, "Expected 3 columns"
        for k,v in df[[columns[0], columns[1]]].values:
            ik2smi[k] = v
        for k,v in df[[columns[0], columns[2]]].values:
            ik2act[k] += [v]
        ik2act = {k: int(np.max(v)) for k,v in ik2act.items()}
        R = []
        for k,v in ik2act.items():
            R += [[k, ik2smi[k], v]]
        return pd.DataFrame(R, columns=columns)

    # Disambiguate data
    all_datasets = {k: disambiguate_data(v) for k,v in all_datasets.items()}
    summary_raw_tasks = []

    print("Printing tasks before last filtering...")
    for dt in sorted(all_datasets):
        l = len(all_datasets[dt])
        columns = list(all_datasets[dt].columns)
        if l != 0:
            ratio = round(sum(all_datasets[dt][columns[2]].tolist()) / l, 3)
        else:
            ratio = 0
        print("{0} -- {1} -- {2}".format(dt, str(l), str(ratio)))

    for k,v in all_datasets.items():
        if v.shape[0] < MIN_SIZE_ANY_TASK:
            continue
        columns = list(v.columns)
        assert len(columns) == 3, "Expected 3 columns"
        n = v[columns[2]].sum()
        if n < MIN_POSITIVES:
            continue
        # if n / len(v) < 0.5:
        file_name = os.path.join(tasks_dir, "{0}_SINGLE_TARGET_{1}.csv".format(k, label))
        print("Saving data in {0}".format(file_name))
        v.to_csv(file_name, index=False)
        summary_raw_tasks.append([k, f'SINGLE TARGET - {label}', len(v), n])

    # Store summary file
    summary_raw_tasks = pd.DataFrame(summary_raw_tasks, columns=["task_id", "target_type", "num_molecules", "num_positives"])
    summary_raw_tasks.to_csv(os.path.join(data_dir, pathogen_code, f"013b_raw_tasks_SP_summary_{label}.csv"), index=False)