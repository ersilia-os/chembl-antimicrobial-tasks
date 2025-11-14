# README for scripts

The `scripts` folder contains scripts that are supposed to run sequentially.

## Step 00. Fetching ChEMBL data

Before running `00_export_chembl_activities.py`, create a **local copy of the ChEMBL DB** following the instructions in `docs/install_ChEMBL.md`. 

Currently, this is set to *ChEMBL_36* (latest version), but if a new release appears simply download the new one and change the DB name in `src/default.py`.

By running `00_export_chembl_activities.py`, a folder named `chembl_activities` will be created inside `config`, containing unmodified activity data extracted directly from ChEMBL tables:

- `activities.csv`
- `assays.csv`
- `assay_parameters.csv`
- `activity_stds_lookup.csv`
- `bioassay_ontology.csv`
- `compound_structures.csv`
- `docs.csv`
- `molecule_dictionary.csv`
- `target_dictionary.csv`

For further information about the content of these tables, please check the official [ChEMBL schema documentation](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.txt).

In addition to raw exports, the script produces three curated files derived from the `ACTIVITIES` table and a simplified version of `ASSAYS`:

- `assay_descriptions.csv`: Table mapping assay IDs to their corresponding text descriptions and ChEMBL IDs.
- `standard_text.csv`: Frequency table of the text stored in the `standard_text_value` column (e.g., 'Compound metabolized')
- `activity_comments.csv`: Frequency table of the comments (text) stored in the `activity_comment` column (e.g., 'active'). 
- `activity_std_units.csv`: Frequency table of the column pairs `standard_type` & `standard_units` (e.g., 'Potency & nM')
- `standard_units.csv`: Frequency table of the standard_units field from the activities.csv file.

After generating the frequency tables (`activity_comments.csv`, `standard_text.csv`, `activity_std_units.csv`, and `standard_units.csv`), a **manual curation process** was performed to assign an activity label to each unique entry. Each item was reviewed and flagged with:

- 1 → indicates the item corresponds to an active compound
- -1 → indicates the item corresponds to an inactive compound
- 0 → indicates the item is inconclusive or ambiguous

## Step 01. Processing compounds

This script calculates the molecular weight (MW) of each compound based on its SMILES (Simplified Molecular Input Line Entry System) representation using RDKit. Running `01_get_compound_info.py` creates the `config/chembl_processed` folder, containing a newly generated file:

- `01_get_compound_info.py`: Table with molregno, chembl_id, molecule_type, canonical_smiles, and calculated MW.

⏳ ETA: ~90 minutes.


## Step 02. Merging activities, assays, compounds & targets

Run `02_merge_activity_tables.py` to merge activity, assay, target, and compound information into a single table, producing the file `activities_all_raw.csv` in `config/chembl_processed`.

⏳ ETA: ~20 minutes.

---
---
---

In this script, ChEMBL tables are automatically curated using an LLM to process activity comments, among others. A series of modified files will be saved in `config/chembl_processed`, each having a column that represents the outcome of the assay/activity (1: Active, -1: Not Active, 0: Unresolved; or Direction: higher value == higher activity --> 1, lower value == higher activity --> -1, inconclusive --> 0). In brief, 4 additional files are generated in  `config/chembl_processed`:

- `activity_comments.csv`: Activity comments are classified as active (1), inactive (-1) or inconclusive (0). New column: `activity_classified`.
- `activity_stds_lookup.csv`: Activity standards (e.g., IC50) are classified with the right direction: the lower the value the higher the activity (-1), the higher the value the higher the activity (1) or inconclusive (0). New column: `activity_direction`.
- `activity_std_units_with_3_assay_descriptions.csv`: Pairs of standard type and unit (e.g., IC50-nM) are provided with 3 randomly sampled descriptions from associated assays, which will be used in the following step.
- `activity_std_units.csv`: Pairs of activity type and unit (e.g., IC50-nM + 3 assay descriptions) are classified with the right direction: the lower the more active (-1), the higher the more active (1) or inconclusive (0). New column: `activity direction`.
- `standard_text.csv`: Text comments are manually processed and classified as active/inactive based on the identification of specific words (). New column: `standard_text_classification`.

To assess the global performance of the LLM in these tasks, we compare the final outcome against a set of manually curated endpoints from ...

# Step 
Uses the LLM processed data to assign an outcome to each activity in ChEMBL (which are summarised in the Activities table).


---
For later on (around step 07):

For this step, it is necessary to download and install Gemma-3 (4b) using Ollama (v0.12.3, check [Ollama's documentation](https://ollama.com/library/gemma3)). The availability of GPUs significantly speeds up this stage. 