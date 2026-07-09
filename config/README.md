# Config files

This folder contains files that have been manually created or curated. You can edit them to tweak or improve the pipeline, but they constitute a base file for further processing.

## Units
We have two files related to unit processing:
- `UnitStringValidations.csv`: maps ChEMBL original units to UCUM-compliant formats. This file was created externally using the [UCUM-LHC](https://ucum.org/) Online Validator and Converter under the _Validate all unit expressions in a CSV file_ section, uploading the file `standard_units.csv` (produced in Step 00 and found in `data/chembl_activities`) and indicating `standard_units` as the name of the column containing the expressions to be validated.
- `ucum_manual.csv`: manual curated effort accounting for more than 290 units (found in `config/ucum_manual.csv`), not only mapping ChEMBL units to valid UCUM formats but also converting units from the same kind to a reference one.

## Synonyms
A single file `synonyms.csv` where the synonyms of the same assay_type (MIC90 & MIC>=90) are curated for simplification.

## Text comments
`standard_text_manual_curation.csv` assigns +1, -1 or 0 to each text value (for example, Not Active) according to active (1), inactive (-1) or unknown (0). Created by manual curation of `data/chembl_activities/standard_text.csv`

## Activity direction
`activity_std_units_manual_curation.csv` assigns the direction of each assay coupled to its unit of measurement (standardised by UCUM). Created by manual processing of `output/stage1/01_activity_std_units_converted.csv`

## Activity comments
`activity_comments_manual_curation.csv` assigns +1,-1 or 0 to each comment (proportion of each comment in the document is shown) according to active (1), inactive (-1) or unknown (0). Created by manual curation of `data/chembl_activities/activity_comments.csv`

## Stage-4 binarisation cutoffs (step 22)
Three files control the active/inactive cutoff applied in `scripts/22_binarise_datasets.py`:
- `category_cutoffs.csv`: the **global default** cutoff per unit, applied to every pathogen (`umol.L-1 → 10`, `% → 50`, `mm → 20`).
- `category_cutoff_candidates.csv`: the **low/middle/high candidate values** per unit (`umol.L-1 5;10;20`, `% 50;75;90`, `mm 10;20;30`) — the only values a per-pathogen override may select. Stringency is direction-dependent: for `umol.L-1` (direction −1) a **lower** value is the more stringent cutoff; for `%`/`mm` (direction +1) a **higher** value is.
- `category_cutoffs_overrides.csv`: **per-pathogen overrides** as `(pathogen_code, unit, tier)`, `tier ∈ {low, middle, high}`. It **ships empty**, so by default every pathogen uses the global default. A per-pathogen cutoff is a threshold decision — fill this in only after reviewing the cutoff-sensitivity diagnostic (`output/stage4/general_plots/27_cutoff_sensitivity.png` / `.csv`, and each pathogen's `22_cutoff_sensitivity.csv`) and with sign-off.

