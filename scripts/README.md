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

<!-- - `activity_std_units.csv`.

    Labels refer to the direction of biological activity:
  - **-1** → lower value = more active (e.g. IC50)
  - **1** → higher value = more active (e.g. %Inhibition)
  - **0** → unclear -->

  - **1** → active
  - **-1** → inactive
  - **0** → inconclusive or unclear

This manual curation allows downstream scripts to automatically use a standardized direction of biological activity for text entries, ensuring consistency across diverse assays and readouts. Both manually curated files are located in `config/manual_curated` (named `activity_comments_manual_curation.csv` and `standard_text_manual_curation.csv`, new column name: `manual_curation`). The user is encouraged to extend or complete these files as needed.

⏳ ETA: ~10 minutes.

## Step 01. Processing compounds

This script calculates the molecular weight (MW) of each compound based on its SMILES (Simplified Molecular Input Line Entry System) representation using RDKit. Running `01_get_compound_info.py` creates the `config/chembl_processed` folder, containing a newly generated file named `compound_info.csv` (table with `molregno`, `chembl_id`, `molecule_type`, `canonical_smiles`, and `calculated MW`).

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

11. **Creating and manually annotating curated [activity type - unit] pairs**: Creating a frequency table of the previously curated column pairs `activity_type` & `unit` (e.g., 'POTENCY & umol·L-1') named `activity_std_units_curated.csv` and located in `config/chembl_processed`. Whenever this last step is done, the user is expected to manually annotate the biological direction of multiple `activity_type` & `unit` pairs, and place the results in a file named `activity_std_units_curated_manual_curation.csv` under the `config/manual_curation` directory (column name: `manual_curation`). This annotation process is conceptually analogous to the ones performed in step 00 to create `activity_comments_manual_curation.csv` and `standard_text_manual_curation.csv`.

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
- `assays_raw.csv`: List of raw assays with metadata (e.g., unit, activity type, compound count).

⏳ ETA: ~4 hours.

## Step 08. Cleaning individual pathogen data

The script `06_clean_pathogen_activities.py` cleans organism-specific activity records for the selected pathogens and summarizes their associated assays. The cleaning steps are enumerated below:



2. **Removing null activities**. Discarding activities with no numerical value nor active or inactive flag in the `activity_comment` nor `standard_text` fields. 

3. **Identifying canonical units**. For each activity type, identifying canonical units i.e. the most occurrying unit. 

Outputs are saved in the folder: `output/<pathogen_code>/`, and include:
- `<pathogen_code>_ChEMBL_cleaned_data.csv`: Cleaned ChEMBL activity records for the selected pathogen (based on fields `target_organism` and `assay_organism`).
- `activity_type_unit_pairs.csv`: List of unique activity type - unit pairs per pathogen, including counts and a canonical unit flag. 
- `assays_cleaned.csv`: List of cleaned assays with metadata (e.g., unit, activity type, compound count).

## Step 07. Clustering assay compounds

The script `07_get_assay_clusters.py` needs to be executed with a conda environment having [bblean](https://github.com/mqcomplab/bblean) installed. 

## Step 08. Assessing compound overlap among assays

## Step 09. Preparing assay data

...

`expert_cutoffs.csv` was manually created with the notebook `expert_cutoffs.ipynb` and can be changed at will. 

## Step 10. Getting compound descriptors

The script `10_get_compound_descriptors.py` needs to be executed with a conda environment having [lazyqsar](hhttps://github.com/ersilia-os/lazy-qsar) installed. 


ETA values were dervied using an ASUS ... 16 CPUs and 32 GB RAM. 


