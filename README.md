# Antimicrobial binary ML tasks from ChEMBL

This repository provides a reproducible pipeline to extract antimicrobial bioactivity data from [ChEMBL](https://www.ebi.ac.uk/chembl/) and convert it into ready-to-use binary classification datasets for machine learning. For each pathogen of interest, the pipeline selects high-quality assays, binarizes quantitative measurements using expert-defined cutoffs, and trains baseline models to evaluate dataset quality, modelability and redundancy. It is aimed at researchers building or benchmarking antimicrobial activity prediction models.

## Setup

To get started, first clone this repository:

```sh
git clone https://github.com/ersilia-os/chembl-antimicrobial-tasks.git
cd chembl-antimicrobial-tasks
```

We recommend creating a Conda environment to run this code. Dependencies are minimal.

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

## Supported pathogens

The pipeline ships with pre-configured settings for 15 pathogens. New pathogens can be added manually by providing a new entry in `config/pathogens.csv` and the corresponding assay allowlist under `config/assays/`.

| Code | Organism |
|------|----------|
| `abaumannii` | *Acinetobacter baumannii* |
| `calbicans` | *Candida albicans* |
| `campylobacter` | *Campylobacter* spp. |
| `ecoli` | *Escherichia coli* |
| `efaecium` | *Enterococcus faecium* |
| `enterobacter` | *Enterobacter* spp. |
| `hpylori` | *Helicobacter pylori* |
| `kpneumoniae` | *Klebsiella pneumoniae* |
| `mtuberculosis` | *Mycobacterium tuberculosis* |
| `ngonorrhoeae` | *Neisseria gonorrhoeae* |
| `paeruginosa` | *Pseudomonas aeruginosa* |
| `pfalciparum` | *Plasmodium falciparum* |
| `saureus` | *Staphylococcus aureus* |
| `smansoni` | *Schistosoma mansoni* |
| `spneumoniae` | *Streptococcus pneumoniae* |

## Pipeline overview

The pipeline has two stages. The first stage processes the raw ChEMBL database into a clean, standardised activity table — this only needs to be run once and takes a few hours. The second stage is pathogen-specific: it filters assays for a given organism, binarizes activities, merges small assays, and trains baseline models to evaluate dataset quality. See [`scripts/README.md`](scripts/README.md) for detailed documentation of every step.

| Stage | Steps | Command |
|-------|-------|---------|
| Process ChEMBL (run once) | 00–06 | `bash scripts/process_ChEMBL.sh` |
| Generate pathogen datasets | 07–21 | `bash scripts/generate_datasets.sh --<code>` |

Runtime for the pathogen stage varies considerably depending on the amount of data available for each organism. For example, to generate datasets for *M. tuberculosis*:
```sh
bash scripts/generate_datasets.sh --mtuberculosis
```

> **Note:** Step 09 requires a GPU and a locally running [ollama](https://ollama.com/) instance for LLM-based assay curation.


## Outputs

After running the pathogen stage, results are written to `output/<pathogen_code>/`. Key files:

| File | Description |
|------|-------------|
| `19_final_datasets.zip` | Selected binary datasets, one CSV per assay (SMILES + binary label) |
| `20_general_datasets.zip` | Pooled binary datasets, one CSV per `(activity_type, unit)` pair across all assays |
| `18_assays_master.csv` | Per-assay annotation table with full pipeline status |
| `20_general_datasets.csv` | Summary table for general organism-level models |
| `21_diagnosis.png` | Multi-panel diagnostic figure |
| `21_summary_table.xlsx` | Summary statistics table |
| `13_individual_LM.csv` | Performance metrics for individual-assay baseline models |
| `15_merged_LM.csv` | Performance metrics for merged-assay baseline models |

## Repository architecture

| Directory | Contents |
|-----------|----------|
| `config/` | Configuration files pre-set by the Ersilia team (cutoffs, pathogen list, manual curation files). Can be edited to customise pipeline behaviour. |
| `data/` | Raw ChEMBL tables and processed intermediate files from Stage 1. |
| `output/` | Per-pathogen results from Stage 2. |
| `scripts/` | All pipeline scripts and the two driver shell scripts. |
| `src/` | Shared utility code used across scripts. |
| `notebooks/` | Exploratory and analysis notebooks. |
| `docs/` | Additional documentation. |


## About the Ersilia Open Source Initiative

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more about us, please visit our [GitBook Documentation](https://ersilia.gitbook.io) and our [GitHub profile](https://github.com/ersilia-os/).

