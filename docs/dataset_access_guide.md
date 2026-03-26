# Dataset Access Guide

This guide explains how to access ML-ready datasets from the ChEMBL antimicrobial pipeline **without manual sifting**. The pipeline provides an automated system for dataset access and organization.

## Overview: No Manual Sifting Required

**The pipeline does NOT require manual sifting through `output/datasets`.** Instead, researchers can systematically access datasets through these key files:

- **`17_final_datasets.csv`** - Lists all selected datasets with metadata
- **`18_assays_master.csv`** - Master reference with comprehensive assay information and selection status
- **Organized ZIP archives** - Datasets grouped by type for efficient loading

## Dataset Access Workflow

### Step 1: Query Selected Datasets

After the pipeline completes, find selected datasets using the final selection table:

```python
import pandas as pd

# Load the final dataset selection table
pathogen_code = "mtuberculosis"  # Replace with your pathogen
final_datasets = pd.read_csv(f'output/{pathogen_code}/17_final_datasets.csv')

# View selected datasets only
selected = final_datasets[final_datasets['selected'] == True]
print(f"Number of selected datasets: {len(selected)}")

# Browse by dataset type
individual_datasets = selected[selected['label'] == 'A']  # Individual assay datasets
merged_datasets = selected[selected['label'] == 'B']      # Merged assay datasets

# View dataset metadata
print(selected[['name', 'activity_type', 'unit', 'target_type', 'cutoff', 'auroc', 'cpds', 'positives']].head())
```

### Step 2: Load Dataset Files Using Existing Utilities

Use the utilities in `src/model_utils.py` to load actual dataset files:

```python
import sys
sys.path.append('src')
from model_utils import load_all_gz_csvs_from_zip, load_data_from_zip

# Load all datasets from ZIP archives
datasets_qt = load_all_gz_csvs_from_zip(f'output/{pathogen_code}/datasets/datasets_qt.zip')
datasets_mx = load_all_gz_csvs_from_zip(f'output/{pathogen_code}/datasets/datasets_mx.zip')
datasets_ql = load_all_gz_csvs_from_zip(f'output/{pathogen_code}/datasets/datasets_ql.zip')

print(f"Loaded {len(datasets_qt)} quantitative, {len(datasets_mx)} mixed, {len(datasets_ql)} qualitative datasets")

# Load a specific dataset by filename
filename = "CHEMBL4649948_PERCENTEFFECT_%_qt_50.0.csv.gz"
if filename in datasets_qt:
    dataset = datasets_qt[filename]
    print(f"Dataset shape: {dataset.shape}")
    print(f"Columns: {list(dataset.columns)}")
```

### Step 3: Access Merged Datasets

Merged datasets are stored in the `datasets/M/` directory:

```python
import os
import gzip

# Load merged datasets directly
merged_dir = f'output/{pathogen_code}/datasets/M'
merged_files = [f for f in os.listdir(merged_dir) if f.endswith('.csv.gz')]

merged_datasets = {}
for filename in merged_files:
    filepath = os.path.join(merged_dir, filename)
    with gzip.open(filepath, 'rt') as f:
        merged_datasets[filename] = pd.read_csv(f)

print(f"Loaded {len(merged_datasets)} merged datasets")
print("Merged dataset examples:", list(merged_datasets.keys())[:5])
```

### Step 4: Map Dataset Names to Files

Link dataset names from the selection table to actual files:

```python
from dataset_utils import make_dataset_filename

# For individual datasets, build filename from metadata
def get_individual_dataset_filename(row):
    """Build dataset filename from 17_final_datasets.csv row"""
    if row['dataset_type'] == 'quantitative':
        return make_dataset_filename(row['name'].split('_')[0], row['activity_type'],
                                   row['unit'], 'quantitative', row['cutoff'])
    elif row['dataset_type'] == 'mixed':
        return make_dataset_filename(row['name'].split('_')[0], row['activity_type'],
                                   row['unit'], 'mixed', row['cutoff'])
    else:
        return make_dataset_filename(row['name'].split('_')[0], row['activity_type'],
                                   row['unit'], 'qualitative')

# For merged datasets, filename is directly available
def get_merged_dataset_filename(row):
    """Get merged dataset filename from 17_final_datasets.csv row"""
    return f"{row['name']}.csv.gz"

# Example usage
for idx, row in selected.head().iterrows():
    if row['label'] == 'A':  # Individual dataset
        filename = get_individual_dataset_filename(row)
        print(f"Individual: {row['name']} -> {filename}")
    else:  # Merged dataset
        filename = get_merged_dataset_filename(row)
        print(f"Merged: {row['name']} -> {filename}")
```

## Dataset File Structure

### Dataset Organization
```
output/<pathogen_code>/
├── 17_final_datasets.csv           # Selected datasets with metadata
├── 18_assays_master.csv           # Master reference table
├── datasets/
│   ├── datasets_qt.zip            # Quantitative datasets
│   ├── datasets_ql.zip            # Qualitative datasets
│   ├── datasets_mx.zip            # Mixed datasets
│   └── M/                         # Merged dataset files
└── correlations/                  # Model predictions for analysis
```

### Dataset File Format
All datasets contain these essential columns:
- **`smiles`** - Standardized compound SMILES string
- **`bin`** - Binary activity label (1=active, 0=inactive)
- **`compound_chembl_id`** - ChEMBL compound identifier
- Additional metadata columns vary by dataset type

### Naming Conventions

**Individual Datasets:**
- Format: `{assay_id}_{activity_type}_{unit}_{type}_{cutoff}.csv.gz`
- Example: `CHEMBL4649948_IC50_umol.L-1_qt_10.0.csv.gz`

**Merged Datasets:**
- Format: `{group_name}_{cutoff}.csv.gz`
- Examples: `M_ORG0_10.0.csv.gz`, `M_SP1_5.0.csv.gz`

## Complete Example: Loading Selected Datasets

```python
import pandas as pd
import sys
sys.path.append('src')
from model_utils import load_all_gz_csvs_from_zip

def load_selected_datasets(pathogen_code):
    """Load all selected datasets for a pathogen"""

    # Load selection metadata
    final_datasets = pd.read_csv(f'output/{pathogen_code}/17_final_datasets.csv')
    selected = final_datasets[final_datasets['selected'] == True]

    # Load individual datasets from ZIP archives
    datasets_qt = load_all_gz_csvs_from_zip(f'output/{pathogen_code}/datasets/datasets_qt.zip')
    datasets_mx = load_all_gz_csvs_from_zip(f'output/{pathogen_code}/datasets/datasets_mx.zip')
    datasets_ql = load_all_gz_csvs_from_zip(f'output/{pathogen_code}/datasets/datasets_ql.zip')

    all_datasets = {**datasets_qt, **datasets_mx, **datasets_ql}

    # Load merged datasets
    import os, gzip
    merged_dir = f'output/{pathogen_code}/datasets/M'
    if os.path.exists(merged_dir):
        for filename in os.listdir(merged_dir):
            if filename.endswith('.csv.gz'):
                filepath = os.path.join(merged_dir, filename)
                with gzip.open(filepath, 'rt') as f:
                    all_datasets[filename] = pd.read_csv(f)

    return selected, all_datasets

# Usage
pathogen_code = "mtuberculosis"
selected_meta, all_datasets = load_selected_datasets(pathogen_code)

print(f"Selected datasets: {len(selected_meta)}")
print(f"Available dataset files: {len(all_datasets)}")
print(f"Example dataset shape: {next(iter(all_datasets.values())).shape}")
```

## Advanced Usage

### Filter by Dataset Properties
```python
# Find high-quality large datasets
high_quality = selected[
    (selected['auroc'] > 0.8) &
    (selected['cpds'] > 5000) &
    (selected['positives'] > 100)
]

# Find specific activity types
ic50_datasets = selected[selected['activity_type'] == 'IC50']
organism_datasets = selected[selected['target_type'] == 'ORGANISM']
```

### Dataset Statistics
```python
# Summary statistics
print("Dataset statistics:")
print(f"Total compounds range: {selected['cpds'].min():.0f} - {selected['cpds'].max():.0f}")
print(f"AUROC range: {selected['auroc'].min():.3f} - {selected['auroc'].max():.3f}")
print(f"Activity types: {selected['activity_type'].unique()}")
print(f"Target types: {selected['target_type'].unique()}")
```

This automated system eliminates the need for manual dataset sifting and provides researchers with direct, programmatic access to ML-ready antimicrobial datasets.