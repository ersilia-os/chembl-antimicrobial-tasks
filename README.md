# Antimicrobial binary ML tasks from ChEMBL

Get antimicrobial tasks from ChEMBL framed as binary classifications. This repository is the updated version of [chembl-binary-tasks](https://github.com/ersilia-os/chembl-binary-tasks). This repository is currently WORK IN PROGRESS.

## Setup

To get started, first clone this repository:

```sh
git clone https://github.com/ersilia-os/chembl-antimicrobial-tasks.git
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

Several configuration data files are needed before gathering and binarizing ChEMBL data:

- activities.csv: 
- activity_comments.csv:
- activity_std_units.csv
- activity_std_units_with_3_assay_descriptions.csv:
- activity_stds_lookup.csv:
- all_activities.csv:
- all_molecules.csv:
- assay_descriptions.csv:
- chembl_35_chemreps.txt:
- pathogens.csv:
- standard_text.csv:
- standard_units_conversions.csv:


🟡 Add brief explanation of each file - once data/README.md is ready!

To download such data, you can use Git LFS:




or download the zip file from Zenodo 🟡 

Alternatively, we provide the code to generate such configuration data. To do it, simply execute:

```bash
sh scripts/00_prepare_config.sh
```

This bash script consecutively executes 4 python scripts extensively described in our [documentation](docs/src_info.md).




### Specifying parameters

Please check *Installing ChEMBL*

PCHEMBL_CUTOFFS = [7]
PERCENTAGE_ACTIVITY_CUTOFFS = [90]
PERCENTILES = [1]
MIN_SIZE_ASSAY_TASK = 500
MIN_SIZE_ASSAY_SUBTASK = 99
MIN_SIZE_ANY_TASK = 99
MAX_NUM_INDEPENDENT_ASSAYS = 5
MAX_NUM_ASSAY_SUBTASKS = 3
MIN_POSITIVES = 10
MAX_NUM_INDEPENDENT_UNITS = 5
DATASET_SIZE_LIMIT = 1e6



## Creating datasets

By default, following scripts assume that PostgreSQL is running in the local machine, and that the database user `chembl_user` with password `1234` has read access to the tables of ChEMBL (database name: `chembl_35`). This can be changed in , as well as the output path the data will be stored in, the minimum number of assays to consider and several parameters for binarization. Feel free to manually edit such variables and values as needed.  