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
`activity_std_units_manual_curation.csv` assigns the direction of each assay coupled to its unit of measurement (standardised by UCUM). Created by manual processing of `data/chembl_processed/00b_activity_std_units_converted.csv`

## Activity comments
`activity_comments_manual_curation.csv` assigns +1,-1 or 0 to each comment (proportion of each comment in the document is shown) according to active (1), inactive (-1) or unknown (0). Created by manual curation of `data/chembl_activities/activity_comments.csv`

