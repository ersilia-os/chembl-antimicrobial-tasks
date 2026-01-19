Before running any script, create a **local copy of the ChEMBL Database** following the instructions detailed in `docs/install_ChEMBL.md`. Currently, this is set to *ChEMBL_36* (latest version), but if a new release appears simply download the new version and change the Database name in `src/default.py` (DATABASE_NAME).

# `scripts` directory

The `scripts` folder contains two scripts that are meant to be run sequentially. 

The first one, named `process_chembl.sh` fetches tables from ChEMBL (step 00), organizes compound data (step 01), standardizes compounds following a well defined and reproducible protocol (step 02), merges heterogeneous data from multiple sources (step 03), converts units (step 04), harmonizes multiple fields (step 04) and creates a single, cleaned file for ChEMBL bioactivities including compound, assay and target information (step 05). This script is meant to be executed only once, and should take no longer than 4-5 hours on a desktop machine. The main computational bottleneck is the standardization of all ChEMBL compounds (step 02). Further details on each step are provided in the `README.md` file inside `src`. 

Files to review before processing ChEMBL:

- `config/activity_comments_manual_curation.csv`: Manual annotation of activity comments generated in step 00 and needed in step 05. 
- `config/standard_text_manual_curation.csv`: Manual annotation of activity standard texts generated in step 00 and needed in step 05.
- `config/synonyms.csv`: Activity type mappings to their associated synonyms. Generated before runnning the pipeline and needed in step 05.
- `config/ucum_GT.csv`: Manual annotation of unit conversions. Generated before runnning the pipeline and needed in step 04.
- `config/UnitStringValidations.csv`: ChEMBL's original unit mappings to UCUM-compliant formats. Generated before runnning the pipeline and needed in step 04.

To process ChEMBL, simply run:

```
conda activate camt
bash process_ChEMBL.sh
```
Optionally, the pipeline can generate Extended Connectivity Fingerprints for all ChEMBL compounds, to be used in subsequent analyses. 
```
bash process_ChEMBL.sh --calculate_ecfps
```

Once a processed version of ChEMBL is available locally, the user can create binarized activity datasets for a given pathogen of interest with `create_datasets.sh`. In brief, the script selects ChEMBL bioactivities associated the the specified pathogen (step 07, )

Files to review before generating binarized datasets for a given pathogen:

- 

```
bash generate_datasets.sh
...
```


