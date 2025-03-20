# Antimicrobial binary ML tasks from ChEMBL

Get antimicrobial tasks from ChEMBL framed as binary classifications. This repository is the updated version of [chembl-binary-tasks](https://github.com/ersilia-os/chembl-binary-tasks).

## Setup

To get started, first clone this repository:

```sh
git clone https://github.com/ersilia-os/chembl-antimicrobial-tasks.git
cd chembl-antimicrobial-tasks
```

We recommend creating a Conda environment to run this code:

```sh
conda create -n adda4tb python=3.10
conda activate adda4tb
pip install -r requirements.txt
```

Access to a postgreSQL database server containing the ChEMBL database is required. You may install ChEMBL in your own computer by following these [instructions](install_ChEMBL.md). To check if the postgreSQL service with the ChEMBL database is up and accessible, you can run the following code with your username, password and database name:

```sh
python src/000_test_db.py --username YOUR_USERNAME --password YOUR_PASSWORD --db_name DB_NAME --organism YOUR_ORGANISM
```

If a data table is output from the script, the checks have been successful. 

## Creating datasets

By default, following scripts assume that PostgreSQL is running in the local machine, and that the database user `chembl_user` with password `1234` has read access to the tables of ChEMBL (database name: `chembl_35`). This can be changed in `src/000_default_parameters.py`, as well as the output path the data will be stored in, the minimum number of assays to consider and several parameters for binarization. Feel free to manually edit such variables and values as needed.  