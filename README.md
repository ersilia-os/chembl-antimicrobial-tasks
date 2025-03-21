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

The following scripts assume that PostgreSQL is running locally, with the username, password, and database name configured in the same file. 


## Creating datasets

Bla bla 