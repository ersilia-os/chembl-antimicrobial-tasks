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
conda create -n camt python=3.12 -y
conda activate camt
pip install -r requirements.txt
```

### Data download
Data is stored with the [eosvc](https://github.com/ersilia-os/eosvc) tool. Please install it in your environment and run both commands to download all data. If you only want access to the cleaned outputs, download only this folder to prevent incurring storage and internet traffic costs:
```sh
eosvc download --path data
eosvc download --path output
```

## Pipeline overview

The pipeline is divided into two main stages: (i) processing ChEMBL (steps 00-06) and (ii) creating binarized bioactivity datasets and training ML models for a given pathogen of interest (steps 07-18). Look at the README file under `/scripts` for detailed information.


## Repository architecture

- `config`: config files pre-set by the Ersilia team. They can be manually edited to obtain different pipeline outputs, if requred
- `data`: all ChEMBL raw and processed data.
- `output`: per-pathogen results


## About the Ersilia Open Source Initiative 🌍🤝

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more about us, please visit our [GitBook Documentation](https://ersilia.gitbook.io) and our [GitHub profile](https://github.com/ersilia-os/).

