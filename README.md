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
python src/00_test_db.py --username YOUR_USERNAME --password YOUR_PASSWORD --db_name chembl_35 --organism YOUR_ORGANISM
```

The code is ... to provide ...