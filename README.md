# Antimicrobial binary ML tasks from ChEMBL

Get antimicrobial tasks from ChEMBL framed as binary classifications. This repository is the updated version of [chembl-binary-tasks](https://github.com/ersilia-os/chembl-binary-tasks). This repository is currently WORK IN PROGRESS.

## Setup

To get started, first clone this repository, avoiding large LFS-stored files:

```sh
GIT_LFS_SKIP_SMUDGE=1 git clone https://github.com/ersilia-os/chembl-antimicrobial-tasks.git
cd chembl-antimicrobial-tasks
```

We recommend creating a Conda environment to run this code. Dependencies are minimal.

```sh
conda create -n camt python=3.10
conda activate camt
pip install -r requirements.txt
```

### Installing ChEMBL

Access to a postgreSQL database server containing the ChEMBL database is required. You may install ChEMBL in your own computer by following these [instructions](docs/install_ChEMBL.md). To check if the postgreSQL service with the ChEMBL database is up and accessible, you can run the following code with your username, password and database name:

```sh
sudo service postgresql start
PGPASSWORD=YOUR_PASSWORD psql -h localhost -p 5432 -U YOUR_USERNAME -d YOUR_DB_NAME -c "\dt"
```
If a List of Relations is displayed, checks have been successfull! Make sure to adapt the variables CHEMBL_USR, CHEMBL_PWD and DATABASE_NAME in `src/default_parameters.py` with your username, password and database name, respectively.

### Downloading configuration data

Several configuration data files are needed before gathering and binarizing ChEMBL data, all of them documented [here](data/README.md). You can pull such data using Git LFS:

```bash
git lfs pull --include="data"
```

Alternatively, we provide the code to generate these data. To do it, simply execute:

```sh
bash scripts/00_prepare_config.sh
```

This bash script consecutively executes 4 Python files extensively described in our [documentation](docs/src_info.md).


### Specifying parameters

We set many parameters to process and binarize ChEMBL bioactivity data, all of which are defined in `src/default_parameters.py`. 

The following scripts assume that PostgreSQL is running locally, with the username, password, and database name configured in the same file. Parameters for binarization are also specified herein. 


## Creating datasets

The primary goal of this repository is to automatically get microbial tasks from ChEMBL framed as a binary classification. To do it, for each pathogen of interest, execute:

```bash
bash scripts/01_fetch_pathogen_data_from_chembl.sh --pathogen_code YOUR_PATHOGEN_CODE --output_dir YOUR_OUTPUT_DIR
```

Note that available pathogen codes are listed in `data/pathogens.csv`. The bash script consecutevely executes 6 Python scripts briefly described as follows:

- `011_pathogen_getter.py`: Retrieves pathogen-related bioactivity data from the ChEMBL database, processes and filters the data, and saves it into structured CSV files for further analysis.
- `012_clean_fetched_pathogen_data.py`:  Reads raw data, applies unit conversions, standardizes activity values, filters relevant information, computes pChEMBL values, and outputs a cleaned dataset in CSV format for further analysis.
- `013_binarize_fetched_pathogen_data.py`: Processes pathogen assay data and binarizes it into datasets using different criteria for machine learning models. Datasets may correspond to specific assays or targets, pChEMBL values, % of activity or global percentiles (sorted by priority).
- `014_datasets_modelability.py`: Computes molecular fingerprints, trains a Random Forest classifier using stratified cross-validation, and evaluates dataset modelability by calculating AUROC scores for each task.
- `015_select_tasks.py`: Selects 25 modelable tasks based on AUROC scores, positive sample ratios, and overlap filtering.
- `016_wrapup_tasks_and_clean_output_folder.py`: Organizes selected tasks into a new directory and creates 2 summary files.

### Output

Many files will be generated when creating the ChEMBL tasks/datasets. Overall, the most important files are:

- `011_{YOUR_PATHOGEN_CODE}_original.csv`: Compounds extracted from ChEMBL and associated to the pathogen of interest. Includes compound information, bioactivity data, assay details, and related metadata. 
- `011_{YOUR_PATHOGEN_CODE}_cleaned.csv`: A cleaned and processed version of the original dataset.
- `014_modelability.csv`: Modelability for each task. Includes AUROC scores to evaluate how well a binary classification model can be trained. Higher AUROCs indicate higher modelability. Tasks have been enumerated on the basis of the parameters specified in `src/default_parameters.py`. 
- `013_raw_tasks and 016_tasks:` For each task, list of active (1) and inactive (0) compounds. `013_raw_tasks` includes all tasks; `016_tasks` includes only the TOP-25 modelable tasks. 
- **`016_tasks_summary.csv`**: Summary of the TOP-25 modelable tasks, accompanied by aggregated statistics and evaluation metrics. 
- **`016_{YOUR_PATHOGEN_CODE}_summary.csv`**: Summary of the final selected tasks specific to the pathogen of interest. 


## About the Ersilia Open Source Initiative

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more about us, please visit our [GitBook Documentation](https://ersilia.gitbook.io) and our [GitHub profile](https://github.com/ersilia-os/).

