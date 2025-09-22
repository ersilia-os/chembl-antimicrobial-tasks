# 🦠 Antimicrobial binary ML tasks from ChEMBL 💊

Get antimicrobial tasks from ChEMBL framed as binary classifications.

This repository is currently **WORK IN PROGRESS**. ⚠️🚧

## Setup 🛠️

To get started, first clone this repository:

```sh
git clone https://github.com/ersilia-os/chembl-antimicrobial-tasks.git
cd chembl-antimicrobial-tasks
```

We recommend creating a Conda environment to run this code. Dependencies are minimal. 🐍

```sh
conda create -n camt python=3.10
conda activate camt
pip install -r requirements.txt
```

<!--

### Installing ChEMBL 🗃️

Access to a postgreSQL database server containing the ChEMBL database is required. You may install ChEMBL in your own computer by following these [instructions](docs/install_ChEMBL.md). To check if the postgreSQL service with the ChEMBL database is up and accessible, you can run the following code with your username, password and database name:

```sh
sudo service postgresql start
PGPASSWORD=YOUR_PASSWORD psql -h localhost -p 5432 -U YOUR_USERNAME -d YOUR_DB_NAME -c "\dt"
```
✅ If a *List of Relations* is displayed, checks have been successfull! ⚠️ Make sure to adapt the variables CHEMBL_USR, CHEMBL_PWD and DATABASE_NAME in `src/default_parameters.py` with your username, password and database name, respectively.

### Downloading configuration data ⚙️

Several configuration data files are needed before gathering and binarizing ChEMBL data, all of them documented [here](data/README.md). You can pull such data using Git LFS:

```bash
git lfs pull --include="data"
```

Alternatively, we provide the code to generate these data. To do it, simply execute:

```sh
bash scripts/00_prepare_config.sh
```

This bash script consecutively executes 4 Python files extensively described in our [documentation](docs/src_info.md).


### Specifying parameters 🧾

We set many parameters to process and binarize ChEMBL bioactivity data, all of which are defined in `src/default_parameters.py`. 

The following scripts assume that PostgreSQL is running locally, with the username, password, and database name configured in the same file. Parameters for binarization are also specified herein. 


## Creating datasets 🔍

The primary goal of this repository is to automatically get microbial tasks from ChEMBL framed as a binary classification. To do it, for each pathogen of interest, execute:

```bash
bash scripts/01_fetch_pathogen_data_from_chembl.sh --pathogen_code YOUR_PATHOGEN_CODE --output_dir YOUR_OUTPUT_DIR
```

Note that available pathogen codes are listed in `data/pathogens.csv`, which can be edited manually. The bash script consecutevely executes 6 Python scripts briefly described as follows:

- `011_pathogen_getter.py`: Retrieves pathogen-related bioactivity data from the ChEMBL database, processes and filters the data, and saves it into structured CSV files for further analysis. 
- `012_clean_fetched_pathogen_data.py`:  Reads raw data, applies unit conversions, standardizes activity values, filters relevant information, computes pChEMBL values, and outputs a cleaned dataset in CSV format for further analysis.
- `013a_binarize_fetched_pathogen_data_ORG.py`: Processes phenotypic-based pathogen assay data and organizes it into datasets that are binarized using different criteria for machine learning models (e.g. pChEMBL, %inhibition, etc). Datasets may correspond to specific assays or targets (i.e. the organism itself), global pChEMBL values, % of activity or comprehensive percentiles (sorted by priority). Datasets are created with six different strategies:

    1. Compounds grouped by assays: fixed assay ID. If the assay has multiple activity types, it's split into several datasets.
    2. Compounds grouped by targets: fixed target ID and activity type and units. Assays may differ. 
    3. Compounds grouped by pChEMBL: assumes the target ID is fixed (i.e. the organism) and integrates all pChEMBL data.
    4. Compounds grouped by percentage: assumes the target ID is fixed (i.e. the organism) and integrates all percentage data.
    5. Compounds grouped by percentiles: fixed target ID (i.e. the organism) - integrates percentile data taking all units into account.
    6. Compounds grouped by activity labels: assumes the target ID is fixed (i.e. the organism) and integrates data using the corresponding activity flag.

    Datasets are binarized following 4 different approaches:

    1. pChEMBL cut-offs
    2. pChEMBL percentiles
    3. Percentage cut-offs
    4. Percentage percentiles

    Datasets not satifying the requirements specified in `src/default_parameters.py` or having a proportion of positives > 0.5 are discarded and not reported. 

- `013b_binarize_fetched_pathogen_data_SP.py`: Processes single protein-based pathogen assay data (both "Binding" and "Functional", separately) and organizes it into datasets that are binarized using different criteria for machine learning models (e.g. pChEMBL, %inhibition, etc). Datasets may correspond to specific assays or targets (e.g. a given protein), global pChEMBL values against a speficic protein, % of activity or comprehensive percentiles (sorted by priority). For further information on dataset creation and binarization please see the previous point. IMPORTANT: in this step strategies (1) and (2) are analogous to `013a_binarize_fetched_pathogen_data_ORG.py`. However, strategies (3), (4), (5) and (6) have been adapted to report results in a target-centric manner (i.e. targets are no longer full organisms but single proteins). 

- `014_datasets_modelability.py`: Computes molecular fingerprints, trains a Random Forest classifier using stratified cross-validation, and evaluates dataset modelability by calculating AUROC scores for each task (i.e. discriminate active compounds from inactives). Additionally, store a Random Forest classifier for each task, trained with all task data.

- `015_datasets_distinguishability.py`: Analogous to dataset modelability, but negative compounds are randomly sampled from ChEMBL. Additionally, store a Random Forest classifier for each task, trained with all task data.



<!--

- `016_select_tasks.py`: Selects 25 modelable tasks based on AUROC scores, positive sample ratios, and overlap filtering.
- `017_wrapup_tasks_and_clean_output_folder.py`: Organizes selected tasks into a new directory and creates 2 summary files.
-->

<!-- ### Output 📊

Many files will be generated when creating the ChEMBL tasks/datasets. Overall, the most important files are:

- `011_{YOUR_PATHOGEN_CODE}_original.csv`: Compounds extracted from ChEMBL and associated to the pathogen of interest. Includes compound information, bioactivity data, assay details, and related metadata. Each row corresponds to a given bioactivity measurement. 
- `011_{YOUR_PATHOGEN_CODE}_cleaned.csv`: A cleaned and processed version of the original dataset. Includes pChEMBL values, %Inhibition, etc.
- `013a_raw_tasks_ORG_summary.csv`: Raw list of phenotypic-based tasks (datasets) created for the pathogen of interest. 
- `013a_raw_tasks_ORG directory`: For each phenotypic-based task (dataset), list of compounds and associated binarized bioactivities.
- `013b_raw_tasks_SP_summary_B.csv`: Raw list of target-based (binding) tasks (datasets) created for the pathogen of interest.
- `013b_raw_tasks_SP_summary_F.csv`: Raw list of target-based (functional) tasks (datasets) created for the pathogen of interest.
- `013a_raw_tasks_SP directory`: For each target-based (both binding and functional) task (dataset), list of compounds and associated binarized bioactivities.
- `014_modelability.csv`: Modelability for each task. Includes AUROC scores to evaluate how well a binary classification model discriminates actives from inactives. Higher AUROCs indicate higher modelability.
- `014_models_MOD.csv`: For each task, performance of a binary classification model trained and tested on the full task data.
- `014_models_MOD directory`: For each task, joblib file including the binary classification model mentioned in the immediately preceding file. 
- `015_distinguishability.csv`: Distinguishability for each task. Includes AUROC scores to evaluate how well a binary classification model using randomly sampled ChEMBL compounds as inactives discriminates actives from inactives. Higher AUROCs indicate higher distinguishability. 
- `015_models_DIS.csv`: For each task, performance of a binary classification model trained and tested on the full task data (negatives are randomly sampled from ChEMBL compounds).
- `015_models_DIS directory`: For each task, joblib file including the binary classification model mentioned in the immediately preceding file. -->

<!--
- `013_raw_tasks and 016_tasks:` For each task, list of active (1) and inactive (0) compounds. `013_raw_tasks` includes all tasks; `017_tasks` includes only all modelable or dist .......... tasks. 
- **`017_tasks_summary.csv`**: Summary of the TOP-25 ........... modelable tasks, accompanied by aggregated statistics and evaluation metrics. 
- **`017_{YOUR_PATHOGEN_CODE}_summary.csv`**: Summary of the final selected tasks specific to the pathogen of interest. 
-->

<!-- ## TL;DR 🚩

Bla bla --> 

## About the Ersilia Open Source Initiative 🌍🤝

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more about us, please visit our [GitBook Documentation](https://ersilia.gitbook.io) and our [GitHub profile](https://github.com/ersilia-os/).

