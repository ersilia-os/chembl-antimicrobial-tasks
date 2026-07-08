This project has been financed by Project PID2023-148309OA-I00 funded by MICIU/AEI/10.13039/501100011033 and by ERDF, EU.

<img width="300" alt="miciu_cofinanciado" src="https://github.com/user-attachments/assets/2ee2934c-9b52-43a6-91ce-4dc5581b2f03" />


# Antimicrobial binary ML tasks from ChEMBL

![status](https://img.shields.io/badge/status-ready-brightgreen)

This repository provides a reproducible pipeline to extract antimicrobial bioactivity data from [ChEMBL](https://www.ebi.ac.uk/chembl/) and turn it into ready-to-use binary classification datasets for machine learning. For each pathogen it selects high-quality assays, binarizes measurements against expert-defined cutoffs, trains baseline models to assess dataset quality, and pools related datasets to maximise usable chemical space. It is aimed at researchers building or benchmarking antimicrobial activity prediction models.

## Setup

```sh
git clone https://github.com/ersilia-os/chembl-antimicrobial-tasks.git
cd chembl-antimicrobial-tasks
conda create -n camt python=3.12 -y
conda activate camt
pip install -r requirements.txt
```

Data (raw inputs and results) is stored with [eosvc](https://github.com/ersilia-os/eosvc). Install it, then download the raw ChEMBL tables and/or the cleaned outputs — download only `output` to avoid the large raw `data`:

```sh
eosvc download --path data      # raw ChEMBL tables (large)
eosvc download --path output    # pipeline results
```

## Pipeline

The pipeline runs in four stages; see [`scripts/README.md`](scripts/README.md) for the detail of every numbered step. Stage 1 runs once (~4–5 h); Stages 2–4 are per pathogen (Stage 4 also has cross-pathogen aggregation steps). All results are written under `output/stage1/` … `output/stage4/`. Each stage has a runner that loops `config/pathogens.csv` (pass codes to restrict); run in the `camt` env.

| Stage | Steps | What it does | How to run |
|-------|-------|--------------|-----------|
| 1 — Process ChEMBL (once) | 00–06 | clean & standardise the raw activity table | `bash scripts/run_stage1.sh` |
| 2 — Per-pathogen prep | 07–12 | filter assays, LLM assay curation, ECFP + eos1klk features | `bash scripts/run_stage2.sh` |
| 3 — Conservative pooling | 13–20 | light models, assay master table, diagnosis | `bash scripts/run_stage3.sh` *(opt-in)* |
| 4 — Signal-based pooling | 21–28 | curate, binarise, transfer-pool, cross-validate | `bash scripts/run_stage4.sh` |

`bash scripts/run_all.sh` chains Stage 2 → Stage 4 (the default end-to-end; Stage 3 is a separate opt-in branch).

> **Note:** steps **09** (GPU + local [ollama](https://ollama.com/) for LLM assay curation) and **12** (`ersilia` env, cross-pathogen eos1klk) are **not** run by `run_stage2.sh` — run them manually (the script prints how). Step 09 is required before Stages 3/4; step 12 before the Stage-4 embedding figures.

## Outputs

Two families of ML-ready datasets are produced per pathogen — the **conservative** per-assay datasets (Stage 3, under `output/stage3/<pathogen_code>/`) and the **signal-based pooled** datasets (Stage 4, under `output/stage4/<pathogen_code>/`). See [`docs/dataset_access_guide.md`](docs/dataset_access_guide.md) for programmatic access, and [`scripts/README.md`](scripts/README.md) for per-step outputs.

## Supported pathogens

15 pathogens ship pre-configured. Add more via `config/pathogens.csv` plus an assay allowlist under `config/assays_to_include/` (and, optionally, an exclusion list under `config/assays_to_exclude/`).

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

## About the Ersilia Open Source Initiative

This repository is developed by the [Ersilia Open Source Initiative](https://ersilia.io). Ersilia develops AI/ML tools to support drug discovery research in the Global South. To learn more, visit our [GitBook](https://ersilia.gitbook.io) and [GitHub](https://github.com/ersilia-os/).
