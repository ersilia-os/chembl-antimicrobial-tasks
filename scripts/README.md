# README for scripts

The `scripts` folder contains scripts that are meant to be run sequentially. Essentially, steps can be grouped in data curation and preprocessing (00-09), data binarization (10-XX) ...

## Step 00. Fetching ChEMBL data

Before running `00_export_chembl_activities.py`, create a **local copy of the ChEMBL Database** following the instructions in `docs/install_ChEMBL.md`. 

Currently, this is set to *ChEMBL_36* (latest version), but if a new release appears simply download the new one and change the Database name in `src/default.py` (DATABASE_NAME).

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

In addition to raw exports, the script produces four curated files derived from the `ACTIVITIES` table and a simplified version of `ASSAYS`:

1. `activity_comments.csv`: Frequency table of the comments (text) stored in the `activity_comment` column (e.g., 'active').
2. `activity_std_units.csv`: Frequency table of the column pairs `standard_type` & `standard_units` (e.g., 'Potency & nM')
3. `standard_units.csv`: Frequency table of the `standard_units` column.
4. `standard_text.csv`: Frequency table of the text stored in the `standard_text_value` column (e.g., 'Compound metabolized')
5. `assay_descriptions.csv`: Table mapping assay IDs to their corresponding text descriptions and ChEMBL IDs.

After generating the frequency tables, a **manual curation process** was performed to assign an activity label to the most frequent entries. Items were reviewed and flagged accordingly in `activity_comments.csv` and `standard_text.csv` (labels refer to compound activity):

  - **1** → active
  - **-1** → inactive
  - **0** → inconclusive or unclear

This manual curation allows downstream scripts to automatically use a standardized direction of biological activity for text entries, ensuring consistency across diverse assays and readouts. Both manually curated files are located in `config/manual_curated` (named `activity_comments_manual_curation.csv` and `standard_text_manual_curation.csv`, new column name: `manual_curation`). The user is encouraged to extend or complete these files as needed.

⏳ ETA: ~5 minutes.

## Step 01. Processing compounds

This script merges ChEMBL compound structure information with compound identifiers to generate a single compound table. Specifically, it combines `config/chembl_activities/compound_structures.csv` with `config/chembl_activities/molecule_dictionary.csv` to map each `molregno` to its corresponding `chembl_id`, and saves the result as `compound_info.csv` in `config/chembl_processed/`. Additionally, molecular weight is calculated for all compounds (`canonical_smiles`) and saved as a `MW` column. 

Outputs are saved in the folder: `config/chembl_processed/`, and include:

`compound_info.csv`: Compound table containing `molregno`, structural fields from `compound_structures.csv`, molecular weight (`MW`) and the corresponding `chembl_id`.

⏳ ETA: ~10 minutes.

## Step 02. Saniziting and standardizing compounds

(...)


### Step 03. Calculating compound descriptors

The script `03_get_compound_descriptors.py` calculates ECFPs (radius 3, 2048 bits) for sanitized and standardized compounds (previous step) using RDKit, ommiting failed SMILES and storing results in H5 file format (`output/descriptors.h5`).

⏳ ETA: ~10 minutes.

## Step 04. Merging activities, assays, compounds & targets

Running `04_merge_all.py` to merge `config/chembl_activities/activities.csv`, `config/chembl_activities/assays.csv`, `config/chembl_activities/target_dictionary.csv`, and `config/chembl_processed/compound_info.csv` into a single table. This script will produce the file `activities_all_raw.csv` in `config/chembl_processed` with a predefined set of columns.

⏳ ETA: ~10 minutes.


## Step 05. Unit harmonization and conversion

The script `05_prepare_conversions.py` generates `unit_conversion.csv` inside `config/chembl_processed/`. This file standardizes measurement units found in ChEMBL activities by:

1. Mapping original ChEMBL unit strings (standard_units) to validated UCUM units (e.g., uM to umol.L-1).
2. Assigning a final standard unit per activity type used in downstream tasks (e.g., IC50 → umol.L-1)
3. Defining a conversion formula (when necessary) to adjust numeric values accordingly (e.g., nmol.L-1 → value/1000 umol.L-1)

Before running this script, make sure the file `UnitStringValidations.csv` mapping ChEMBL's original units to UCUM-compliant formats is available in `config/chembl_processed/`. This file was created externally using the [UCUM-LHC](https://ucum.org/) Online Validator and Converter under the _Validate all unit expressions in a CSV file_ section, uploading the file `standard_units.csv` (produced in Step 00 and found in `config/chembl_activities`) and indicating `standard_units` as the name of the column containing the expressions to be validated. These string validations are merged with a manual curation effort accounting for more than 290 units (found in `config/manual_curation/ucum_GT.csv`), not only mapping ChEMBL units to valid UCUM formats but also converting units from the same kind to a reference one. The file `unit_conversion.csv` is created in `config/chembl_processed` and includes all the information mentioned above. 

⏳ ETA: ~0 minutes.

## Step 06. Cleaning activities table

The script `06_preprocess_activity_data.py` produces a curated and standardized version of the activity data, saved as `activities_preprocessed.csv` in `config/chembl_processed`. The file contains a final cleaned and normalized dataset with all compound-assay-target activity records. This step performs several cleaning and harmonization subtasks:

1. **Filtering invalid entries**. Removing activities with missing canonical_smiles (226k). No other entries are removed during this cleaning process. 

2. **Flagging activity comments**. Loading manual annotations from `activity_comments_manual_curation.csv` and flagging each activity comment as active (1), inactive (-1) or unknown (0).

3. **Flagging standard text comments**. Loading manual annotations from `standard_text_manual_curation.csv` and flagging each standard text comment as active (1), inactive (-1) or unknown (0).

4. **Harmonizing and converting units**. Using `unit_conversion.csv` (from Step 05) to normalize unit strings and convert raw values using predefined conversion formulas. It produces `converted_units.csv` (frequencies of new units) and `converted_units_map.csv` (mappings from original to final units), both located in `config/chembl_processed`.

5. **Normalizing activity types**. Normalizing variations in the `standard_type` column (e.g., 'A ctivity', 'Activ ity', 'Activit y', 'Activity', 'activity'). Outputs `harmonized_types_map.csv` for reference.

6. **Normalizing relations**. Mapping relations like ">=", ">>", "<=" to simplified forms (>, <, =).

7. **Recalculating pChEMBL values**. Recalculating pChEMBL values when the unit is umol.L-1 and a numeric value is available.

8. **Replacing Document ID**. Replacing `doc_id` with its corresponding `doc_chembl_id` using the `docs.csv` table from `config/chembl_activities`.

9. **Column renaming and cleanup**. Dropping original fields (e.g., `standard_value`, `standard_units`) and renaming columns for clarity (`value`, `unit`, `relation`, `activity_type`, `activity_comment`, `standard_text`, etc.)

10. **Converting activity types to their corresponding synonyms**. Mapping activity types to their synonyms as defined in `config/manual_curation/synonyms.csv`. 

11. **Creating and manually annotating curated [activity type - unit] pairs**: Creating a frequency table of the previously curated column pairs `activity_type` & `unit` (e.g., 'POTENCY & umol·L-1') named `activity_std_units_curated.csv` and located in `config/chembl_processed`. Whenever this last step is done, the user is expected to manually annotate the biological direction of multiple `activity_type` & `unit` pairs (see examples below), and place the results in a file named `activity_std_units_curated_manual_curation.csv` under the `config/manual_curation` directory (column name: `manual_curation`). In brief, the curated label refers to the direction in which biological activity increases:

  - **-1** → lower value = more active (e.g. IC50)
  - **1** → higher value = more active (e.g. %INHIBITION)
  - **0** → unclear or inconclusive

    Additionally, the script generates a summary table (`activity_std_units_counts_unit_comment.csv`, saved in `config/chembl_processed`) reporting, for each `activity_type`, the number of records with:
    (i) unit + comment, (ii) no unit + comment, (iii) unit + no comment, and (iv) no unit + no comment. A comment is defined as present when either `activity_comment` or `standard_text` is non-zero The table is sorted by total record count per activity type and is intended for quick inspection of quantitative vs qualitative coverage.

⏳ ETA: ~15 minutes.

## Step 07. Splitting data by pathogen

The script `07_get_pathogen_assays.py` filters the full preprocessed ChEMBL dataset (`activities_preprocessed.csv`) to extract organism-specific subsets and summarizes their associated assays.

Currently, the following **list of pathogens** is processed:

```bash
["Acinetobacter baumannii", "Candida albicans", "Campylobacter", "Escherichia coli", "Enterococcus faecium", "Enterobacter", "Helicobacter pylori", "Klebsiella pneumoniae", "Mycobacterium tuberculosis", "Neisseria gonorrhoeae", "Pseudomonas aeruginosa", "Plasmodium falciparum", "Staphylococcus aureus", "Schistosoma mansoni", "Streptococcus pneumoniae"]
 ```

Pathogen codes are automatically generated by taking the first letter of the genus and appending the species name, all lowercase.
For example:

- *Escherichia coli* → ecoli
- *Staphylococcus aureus* → saureus
- *Mycobacterium tuberculosis* → mtuberculosis

Outputs are saved in the folder: `output/<pathogen_code>/`, and include:
- `<pathogen_code>_ChEMBL_raw_data.csv`: All ChEMBL activity records for the selected pathogen (based on fields `target_organism` and `assay_organism`).
- `target_organism_counts.csv`: Frequency of target organisms found in the data.
- `compound_counts.csv`: Frequency of compounds found in the data.
- `assays_raw.csv`: List of raw assays with metadata (e.g., unit, activity type, compound count). Each [`assay_id`, `activity_type`, `unit`] item is treated independently. 

⏳ ETA: ~4 hours. [REVISE]

## Step 08. Cleaning individual pathogen data

The script `08_clean_pathogen_activities.py` cleans organism-specific activity records for the selected pathogens and summarizes their associated assays. The cleaning steps are enumerated below:

1. **Removing null activities**. Discarding activities with no numerical value or active/inactive flag in the `activity_comment` nor `standard_text` fields. 

2. **Identifying activity directions**. For each activity type and unit pair, identifying the direction of the biological activity. 

3. **Removing unmodelable activities**. Keeping only those activities with [-1, +1] direction or active or inactive flag in the `activity_comment` nor `standard_text` fields.

4. **Identifying canonical units**. For each activity type, identifying canonical units i.e. the most occurrying unit for that specific pathogen. 

Outputs are saved in the folder: `output/<pathogen_code>/`, and include:
- `<pathogen_code>_ChEMBL_cleaned_data.csv`: Cleaned ChEMBL activity records for the selected pathogen (based on fields `target_organism` and `assay_organism`).
- `activity_type_unit_pairs.csv`: List of unique activity type - unit pairs per pathogen, including counts, canonical unit flags, and assigned biological direction (when available). 
- `assays_cleaned.csv`: List of cleaned assays with metadata (e.g., unit, activity type, compound count).

## Step 09. Calculating assay clusters

The script `09_calculate_assay_clusters.py` needs to be executed with a conda environment having [bblean](https://github.com/mqcomplab/bblean) installed. For each individual assay ([`assay_id`, `activity_type`, `unit`] item), unique compounds are clustered based on their ECFP4 fingerprints (2048 bits) using the BitBirch algorithm. The number of clusters is computed at three distinct Tanimoto Coefficient cut-offs: 0.3, 0.6, and 0.85, providing a measure of chemical diversity within each assay.

## Step 10. Assessing compound overlap among assays

The script `10_get_assay_overlap.py` computes pairwise compound overlap between assays within each pathogen-specific cleaned dataset. Each assay is treated as an independent [`assay_id`, `activity_type`, `unit`] item. Only assays with 50 or more compounds are considered. For each assay pair, the number of shared compounds and a normalized overlap ratio are calculated.

Outputs are saved in the folder: `output/<pathogen_code>/`, and include:
- assays_overlap.csv: Pairwise assay overlap table, including compound counts, shared compounds, and overlap ratios.

⏳ ETA: [REVISE]

## Step 11. Curating assay parameters

The script `11_curate_assay_parameters.py` extracts and standardizes additional assay-level biological context for each pathogen, using a local LLM (via `ollama`) to parse existing ChEMBL assay annotations and associated publication metadata. For each pathogen, the script iterates over all cleaned assays (`assays_cleaned.csv`) and compiles a text block containing assay fields from ChEMBL (`config/chembl_activities/assays.csv`), document metadata (`config/chembl_activities/docs.csv`), and summary statistics from Step 08 (e.g., compound counts, activity type, unit, direction). This block is passed to an information-extraction prompt, requesting a strict JSON output with the following curated fields: [`organism`, `target_type_curated`, `strain`, `atcc_id`, `mutations`, `known_drug_resistances`, `media`]. One JSON file is generated per [`assay_id`, `activity_type`, `unit`] item and saved under `output/<pathogen_code>/assay_parameters/`. After processing all assays, the folder is compressed into `output/<pathogen_code>/assay_parameters.zip` and the temporary directory is removed.

Outputs are saved in the folder: `output/<pathogen_code>/`, and include:
- `assay_parameters.zip`: ZIP archive containing one *_parameters.json file per assay item with curated assay parameters.

⏳ ETA: [REVISE]


## Step 12: Preparing assay data

The script `12_prepare_assay_data.py` prepares assay-level datasets for each pathogen by combining cleaned activity data (Step 08), curated assay parameters (Step 11), and expert-defined cutoffs. For each [`assay_id`, `activity_type`, `unit`] item in `assays_cleaned.csv`, the script generates binarized quantitative and/or qualitative datasets at the compound level. Main operations include:

1. **Loading curated target type**. Reading `target_type_curated` for each assay item from `output/<pathogen_code>/assay_parameters.zip`.

2. **Loading expert cutoffs**. Reading expert thresholds from `config/manual_curation/expert_cutoffs.csv` (keyed by [`activity_type`, `unit`, `target_type_curated`, `pathogen_code`]).
⚠️ Important: `expert_cutoffs.csv` is manually created using the notebook `expert_cutoffs.ipynb` and can be modified or updated at will.

3. **Quantitative dataset preparation** (when possible). For assays with a valid direction and quantitative values:

- Adjusting censored relations (</>) according to the activity direction.
- Selecting a single measurement per compound (most active value).
- Binarizing values using the expert cutoff (if available).

4. **Qualitative dataset preparation** (when possible). Building compound-level binary labels from `activity_comment` and `standard_text`, enforcing consistency and removing conflicts.

5. **Dataset typing and summary**. Assigning each assay a `dataset_type` (quantitative, qualitative, mixed, or none) and storing summary statistics and label distributions (number of positives, ratio, etc).

Outputs are saved in the folder: **output/<pathogen_code>/**, and include:
- `datasets/`: per-assay datasets saved as compressed CSV files. 
- `*_qt.csv.gz`: binarized quantitative dataset (only if an expert cutoff is available).
- `*_ql.csv.gz`: binarized qualitative dataset (only if qualitative labels exist).
- `assays_data.csv`: assay-level summary table including dataset type, expert cutoff, relation counts, basic activity statistics, and label ratios.

⏳ ETA: [REVISE]

## Step 13. Merging assay tables

The script `13_merge_assay_tables.py` merges all assay-level tables generated in previous steps into a single master table per pathogen. Assays are treated as independent [`assay_id`, `activity_type`, `unit`] items. For each pathogen, the script loads: `assays_cleaned.csv` (Step 08), `assays_clusters.csv` (Step 09), `assay_parameters.zip` (Step 11) and `assays_data.csv` (Step 12). It first checks that each table is unique on the keys [`assay_id`, `activity_type`, `unit`]. Then, curated assay parameters are extracted from `assay_parameters.zip` and appended to `assays_cleaned.csv`. Finally, the cluster metrics and dataset summaries are merged into the same table.

Outputs are saved in the folder: **output/<pathogen_code>/**, and include:
- `assays_master.csv`: Master assay table combining cleaned assay metadata, curated parameters, clustering statistics, and dataset summaries.

⏳ ETA: ~0 minutes.






ETA values were dervied using an ASUS ... 16 CPUs and 32 GB RAM. 


