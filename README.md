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

The primary purpose of this repository is to automatically get microbial tasks from ChEMBL framed as a binary classification. To do it, for each pathogen of interest, execute:

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

