### Detailed information for all python scripts in src/

- `001_units_resolver_using_statistics.py`: This script uses activity annotations from `standard_text.csv` and `activity_comments.csv` to label entries in `activities.csv`, statistically compares active vs. inactive values to infer activity direction per (`standard_type`, `standard_units`), and writes these directions to `activity_std_units.csv`.

    > Inputs: `standard_text.csv`, `activity_comments.csv`, `activities.csv`

    > Outputs: Updated `activities.csv` (with `activity_flag` column), created `activity_std_units.csv` (with `activity_direction` column)


- `002_activity_binarizer_and_directions_with_llm.py`: This script uses rule-based logic and a local LLM (`LLM_BIN_FILENAME`) to classify entries in `activity_comments.csv` as active (1), inactive (-1), or inconclusive (0), and determines the directionality of activity values for each (`standard_type`, `standard_units`) pair in `activity_stds_lookup.csv` and `activity_std_units.csv`, optionally leveraging up to three example assay descriptions per unit to improve classification.

    > Inputs: `activity_comments.csv`, `activity_stds_lookup.csv`, `activity_std_units.csv`, `assay_descriptions.csv`

    > Outputs: Updated `activity_comments.csv` (with `activity_classified` column), updated `activity_stds_lookup.csv` and `activity_std_units.csv` (with `activity_direction` column), created `activity_std_units_with_3_assay_descriptions.csv`

- `003_chembl_chemical_space.py`: This script processes chemical structures from `chembl_35_chemreps.txt` by computing the number of heavy atoms and the InChIKey for each molecule using RDKit, and writes the results to a new CSV file.

    > Inputs: `chembl_35_chemreps.txt`

    > Outputs: Created `all_molecules.csv` (with columns: `chembl_id`, `smiles`, `num_heavy_atoms`, `inchikey`)


- `004_clean_all_activities.py`: This script consolidates and cleans activity data from multiple sources by assigning activity and direction flags, normalizing standard relation values, enriching entries with direction confidence, and generating a curated activity table along with a summary of standard unit usage and conversions.

    > Inputs: `activities.csv`, `activity_comments.csv`, `standard_text.csv`, `activity_stds_lookup.csv`, `activity_std_units.csv`, `standard_units_conversions.csv` (optional)

    > Outputs: Created `all_activities.csv` (cleaned and enriched activity data), created or overwritten `standard_units_conversions.csv` (with usage count, conversion formula, and new unit)