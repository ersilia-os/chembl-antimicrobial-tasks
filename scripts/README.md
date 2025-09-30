# README for scripts

The `scripts` folder contains scripts that are supposed to run sequentially.

# Step 000

Before running `000_export_chembl_activities.py`, create a local copy of the ChEMBL DB following the instructions in `docs/install_ChEMBL.md`. Currently, this is set to ChEMBL_36 (latest version), but if a new release appears simply download the new one and change the DB name in `default.py`.

By running `000_export_chembl_activities.py`, a folder named `chembl_activities` will be created inside `config`, containing unmodified activity data extracted directly from ChEMBL tables. For further information about the content of these tables (`ACTIVITIES`, `ACTIVITY_PROPERTIES`, `ACTIVITY_STDS_LOOKUP`, `ACTIVITY_SUPP`, `ACTIVITY_SUPP_MAP`, `ACTIVITY_SMID`, `ACTION_TYPE`, `ASSAYS`), please check the official [ChEMBL schema documentation](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.txt).

In addition to raw exports, the script produces three curated files derived from the `ACTIVITIES` table and a simplified version of `ASSAYS`:

- `activity_comments.csv`: Frequency table of the comments (text) stored in the `activity_comment` column (e.g., 'active'). 
- `activity_std_units.csv`: Frequency table of the column pairs `standard_type` & `standard_units` (e.g., 'Potency & nM')
- `standard_text.csv`: Frequency table of the text stored in the `standard_text_value` column (e.g., 'Compound metabolized')
- `assay_descriptions.csv`: Table mapping assay IDs to their corresponding text descriptions and ChEMBL IDs.

# Step 001   --> to discuss

In this script, ChEMBL tables are automatically curated using an LLM to process activity comments, among others. A series of modified files will be saved in `config/llm_processed`, each having a column that represents the outcome of the assay/activity (1: Active, -1: Not Active, 0: Unresolved; or Direction: higher value == higher activity --> 1, lower value == higher activity --> -1, inconclusive --> 0). For this step, it is necessary to download and install Gemma-3 (4b) using Ollama (v0.12.3, check [Ollama's documentation](https://ollama.com/library/gemma3)). The availability of GPUs significantly speeds up this stage. In brief, 4 additional files are generated in  `config/llm_processed`:

- `activity_comments.csv`: Activity comments are classified as active (1), inactive (-1) or inconclusive (0). New column: `activity_classified`.
- `activity_stds_lookup.csv`: Activity standards (e.g., IC50) are classified with the right direction: the lower the value the higher the activity (-1), the higher the value the higher the activity (1) or inconclusive (0). New column: `activity_direction`.
- `activity_std_units_with_3_assay_descriptions.csv`: Pairs of standard type and unit (e.g., IC50-nM) are provided with 3 randomly sampled descriptions from associated assays, which will be used in the following step.
- `activity_std_units.csv`: Pairs of activity type and unit (e.g., IC50-nM + 3 assay descriptions) are classified with the right direction: the lower the more active (-1), the higher the more active (1) or inconclusive (0). New column: `activity direction`.
- `standard_text.csv`: Text comments are manually processed and classified as active/inactive based on the identification of specific words (). New column: `standard_text_classification`.

To assess the global performance of the LLM in these tasks, we compare the final outcome against a set of manually curated endpoints from ...

# Step 002
Uses the LLM processed data to assign an outcome to each activity in ChEMBL (which are summarised in the Activities table).