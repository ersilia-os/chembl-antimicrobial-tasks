

## Step 00. Fetching ChEMBL data


By running `00_export_chembl_activities.py`, a folder named `chembl_activities` will be created inside `data`, containing unmodified activity data extracted directly from ChEMBL tables:

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

This manual curation allows downstream scripts to automatically use a standardized direction of biological activity for text entries, ensuring consistency across diverse assays and readouts. Both manually curated files are located in `config` (named `activity_comments_manual_curation.csv` and `standard_text_manual_curation.csv`, new column name: `manual_curation_activity`). The user is encouraged to extend or complete these files at will.

⏳ ETA: ~5 minutes.

## Step 01. Processing compounds

This script merges ChEMBL compound structure information with compound identifiers to generate a single compound table. Specifically, it combines `data/chembl_activities/compound_structures.csv` with `data/chembl_activities/molecule_dictionary.csv` to map each `molregno` to its corresponding `chembl_id`, and saves the result as `compound_info.csv` in `data/chembl_processed/`. Additionally, molecular weight is calculated for all compounds (`canonical_smiles`) and saved as a `MW` column. 

Outputs are saved in the folder: `data/chembl_processed/`, and include:

`compound_info.csv`: Compound table containing `molregno`, structural fields from `compound_structures.csv`, molecular weight (`MW`) and the corresponding `chembl_id`.

⏳ ETA: ~10 minutes.

## Step 02. Saniziting and standardizing compounds

This step standardizes compound structures to ensure consistent molecular representations in downstream tasks. Starting from `data/chembl_processed/compound_info.csv`, the script canonicalizes SMILES, removes salts and solvents, and extracts the parent molecule using the ChEMBL structure pipeline. Molecular weight is then recalculated based on the standardized structure.

Outputs are saved in the folder: `data/chembl_processed/`, and include:

`compound_info_standardized.csv`: Table containing standardized parent SMILES (`standardized_smiles`) and recalculated molecular weight (`standardized_MW`).

⏳ ETA: ~5 hours.


## Step 03. Merging activities, assays, compounds & targets

Running `03_merge_all.py` to merge `data/chembl_activities/activities.csv`, `data/chembl_activities/assays.csv`, `data/chembl_activities/target_dictionary.csv`, and `data/chembl_processed/compound_info.csv` into a single table. This script will produce the file `activities_all_raw.csv` in `data/chembl_processed` with a predefined set of columns.

⏳ ETA: ~10 minutes.


## Step 04. Unit harmonization and conversion

The script `04_prepare_conversions.py` generates `unit_conversion.csv` inside `data/chembl_processed/`. This file standardizes measurement units found in ChEMBL activities by:

1. Mapping original ChEMBL unit strings (standard_units) to validated UCUM units (e.g., uM to umol.L-1).
2. Assigning a final standard unit per activity type used in downstream tasks (e.g., IC50 → umol.L-1)
3. Defining a conversion formula (when necessary) to adjust numeric values accordingly (e.g., nmol.L-1 → value/1000 umol.L-1)

Before running this script, make sure the file `UnitStringValidations.csv` mapping ChEMBL's original units to UCUM-compliant formats is available in `config/`. This file was created externally using the [UCUM-LHC](https://ucum.org/) Online Validator and Converter under the _Validate all unit expressions in a CSV file_ section, uploading the file `standard_units.csv` (produced in Step 00 and found in `data/chembl_activities`) and indicating `standard_units` as the name of the column containing the expressions to be validated. These string validations are merged with a manual curation effort accounting for more than 290 units (found in `config/ucum_GT.csv`), not only mapping ChEMBL units to valid UCUM formats but also converting units from the same kind to a reference one. The file `unit_conversion.csv` is created in `data/chembl_processed` and includes all the information mentioned above. 

⏳ ETA: ~0 minutes.

## Step 05. Cleaning activities table

The script `05_clean_activities.py` produces a curated and standardized version of the activity data, saved as `activities_preprocessed.csv` in `data/chembl_processed`. The file contains a final cleaned and normalized dataset with all compound-assay-target activity records. This step performs several cleaning and harmonization subtasks:

1. **Filtering invalid entries**. Removing activities with missing canonical_smiles (226k). No other entries are removed during this cleaning process. 

2. **Flagging activity comments**. Loading manual annotations from `config/activity_comments_manual_curation.csv` and flagging each activity comment as active (1), inactive (-1) or unknown (0).

3. **Flagging standard text comments**. Loading manual annotations from `config/standard_text_manual_curation.csv` and flagging each standard text comment as active (1), inactive (-1) or unknown (0).

4. **Harmonizing and converting units**. Using `unit_conversion.csv` (generated in Step 04) to normalize unit strings and convert raw values using predefined conversion formulas. It produces `converted_units.csv` (frequencies of new units) and `converted_units_map.csv` (mappings from original to final units), both located in `data/chembl_processed`.

5. **Normalizing activity types**. Normalizing variations in the `standard_type` column (e.g., 'A ctivity', 'Activ ity', 'Activit y', 'Activity', 'activity'). Outputs `harmonized_types_map.csv` for reference.

6. **Normalizing relations**. Mapping relations like ">=", ">>", "<=" to simplified forms (>, <, =).

7. **Recalculating pChEMBL values**. Recalculating pChEMBL values when the unit is umol.L-1 and a numeric value is available.

8. **Replacing Document ID**. Replacing `doc_id` with its corresponding `doc_chembl_id` using the `docs.csv` table from `data/chembl_activities`.

9. **Column renaming and cleanup**. Dropping original fields (e.g., `standard_value`, `standard_units`) and renaming columns for clarity (`value`, `unit`, `relation`, `activity_type`, `activity_comment`, `standard_text`, etc.)

10. **Converting activity types to their corresponding synonyms**. Mapping activity types to their synonyms as defined in `config/synonyms.csv`. 

11. **Creating and manually annotating curated [activity type - unit] pairs**: Creating a frequency table of the previously curated column pairs `activity_type` & `unit` (e.g., 'POTENCY & umol·L-1') named `activity_std_units_curated.csv` and located in `data/chembl_processed`. Whenever this last step is done, the user is expected to manually annotate the biological direction of multiple `activity_type` & `unit` pairs (see examples below), and place the results in a file named `activity_std_units_curated_manual_curation.csv` under the `config` directory (column name: `manual_curation_direction`). In brief, the curated label refers to the direction in which biological activity increases:

  - **-1** → lower value = more active (e.g. IC50)
  - **1** → higher value = more active (e.g. %INHIBITION)
  - **0** → unclear or inconclusive

    Additionally, the script generates a summary table (`activity_std_units_curated_comments.csv`, saved in `data/chembl_processed`) reporting, for each `activity_type` and `unit` pair, the total number of counts as well as the number of records with a non-zero value (1 for active and -1 for inactive) in `activity_comment` or `standard_text`.

⏳ ETA: ~15 minutes.


### Step 06. Calculating Extended Connectivity Fingerprints (OPTIONAL)

The script `06_calculate_ecfps.py` calculates ECFPs (radius 3, 2048 bits) for sanitized and standardized compounds (step 02) using RDKit, ommiting failed SMILES and storing results in H5 file format (`data/chembl_processed/ChEMBL_ECFPs.h5`).

⏳ ETA: ~15 minutes.

## Step 07. Splitting data by pathogen

The script `07_get_pathogen_assays.py` filters the full preprocessed ChEMBL dataset (`data/chembl_processed/activities_preprocessed.csv`) to extract pathogen-specific subsets and summarizes their associated assays. Supported pathogen codes are listed in `config/pathogens.csv`. Bioactivity data is gathered by searching the column `pathogen` from `config/pathogens.csv` in `target_organism` and `assay_organism` fields from `data/chembl_processed/activities_preprocessed.csv`. The user might include valuable ChEMBL assays that miss this condition within `config/assays/<pathogen_code>.csv` (comma- or newline- separated assay ChEMBL IDs).

Outputs are saved in the folder: `output/<pathogen_code>/`, and include:
- `<pathogen_code>_ChEMBL_raw_data.csv`: All ChEMBL activity records for the selected pathogen (based on fields `target_organism` and `assay_organism`).
- `target_organism_counts.csv`: Frequency of target organisms found in the data.
- `compound_counts.csv`: Frequency of compounds found in the data.
- `assays_raw.csv`: List of raw assays with metadata (e.g., unit, activity type, compound count). Each [`assay_id`, `activity_type`, `unit`] item is treated independently. 

⏳ ETA: ~1 minute.

## Step 08. Cleaning individual pathogen data

The script `08_clean_pathogen_activities.py` cleans organism-specific ChEMBL activity records the selected pathogen and produces summarized assay-level outputs. The main steps are:

1. **Removing non standardized compounds**: Activities associated to compounds that have not been standardized when processing ChEMBL (step 02) are discarded. 

2. **Removing empty activities**: Activities lacking both a numerical value and a non-zero `text_flag` are discarded.

3. **Filtering by supported units**: Only activities with consensus units (from `data/chembl_processed/unit_conversion.csv`) or missing units are retained.

4. **Assigning activity direction**: Biological direction (-1, 0, +1) is assigned per (`activity_type`, `unit`) pair using manually curated annotations from `config/activity_std_units_curated_manual_curation.csv`.

5. **Removing unmodelable activities**: Activities are kept only if they have a clear direction (-1 or +1) or an active/inactive `text_flag`.

6. **Summarizing activity type–unit pairs**: Unique (`activity_type`, `unit`) pairs are counted, annotated with direction and text comments and ranked by frequency.

7. **Aggregating cleaned assay data**: Assays are split by `activity_type` and `unit`, summarized with metadata (compound counts, missing values, direction, text annotations).

Outputs are written to **output/<pathogen_code>/** and include:

- `<pathogen_code>_ChEMBL_cleaned_data.csv.gz`: cleaned activity-level data for the selected pathogen (based on fields `target_organism` and `assay_organism`).
- `activity_type_unit_pairs.csv`: counts and directions for activity type–unit combinations.
- `activity_type_unit_comment.csv`: activity type summaries by unit availability and text annotations.
- `assays_cleaned.csv`: cleaned and filtered assay-level metadata.
  

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






ETA values were dervied using an ASUS ... 16 CPUs and 32 GB RAM. Case study mtb.

Storage space.


