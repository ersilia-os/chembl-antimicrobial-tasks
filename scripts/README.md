# Scripts

Two sequential stages: **Stage 1** processes the full ChEMBL database once; **Stage 2** runs per pathogen to produce ML-ready datasets.

## Setup

Create a local ChEMBL PostgreSQL database following `docs/install_ChEMBL.md`. Database connection defaults are in `src/default.py` (`DATABASE_NAME`, user, password).

---

## Stage 1 — Process ChEMBL (run once, ~4–5 h)

```sh
conda activate camt
bash scripts/process_ChEMBL.sh              # standard
bash scripts/process_ChEMBL.sh --calculate_ecfps  # also compute fingerprints
```

| Step | Script | What it does | Key output |
|------|--------|-------------|------------|
| 00 | `00_export_chembl_activities.py` | Exports raw tables from the ChEMBL PostgreSQL DB and derives frequency tables for activity comments, units, and text values. | `data/chembl_activities/*.csv` |
| 01 | `01_prepare_manual_files.py` | Harmonizes activity type strings, maps units to UCUM, collapses synonyms. **Pauses here** — a curator must assign biological direction to each (activity\_type, unit) pair before continuing (see below). | `data/chembl_processed/01_*.csv` |
| 02 | `02_get_compound_info.py` | Merges compound structures with molecule dictionary; calculates molecular weight. | `data/chembl_processed/02_compound_info.csv` |
| 03 | `03_standardize_compounds.py` | Canonicalizes SMILES, removes salts/solvents, recalculates MW. ⏳ ~3–4 h | `data/chembl_processed/03_compound_info_standardized.csv` |
| 04 | `04_merge_activity_and_compounds.py` | Joins activities, assays, targets, and compounds into a single table. | `data/chembl_processed/04_activities_all_raw.csv` |
| 05 | `05_clean_activities.py` | Filters invalid entries, applies activity comment/text flags, converts units, harmonizes activity types, recalculates pChEMBL. ⏳ ~15 min | `data/chembl_processed/05_activities_preprocessed.csv` |
| 06 | `06_calculate_ecfps.py` *(optional)* | Computes ECFP4 fingerprints (radius 3, 2048 bits) for all standardized compounds. ⏳ ~15 min | `data/chembl_processed/06_chembl_ecfps.h5` |

### Manual curation required between steps 01 and 02

1. Open `data/chembl_processed/01_activity_std_units_converted.csv`
2. Fill in `manual_curation_direction`: `1` = higher value → more active, `-1` = lower value → more active, `0` = unclear
3. Save as `config/activity_std_units_manual_curation.csv`

### Config files used in Stage 1

| File | Step | Purpose |
|------|------|---------|
| `config/ucum_manual.csv` | 01, 05 | Maps ChEMBL units to UCUM; provides value conversion formulas |
| `config/synonyms.csv` | 01, 05 | Collapses activity type variants to canonical names |
| `config/activity_comments_manual_curation.csv` | 05 | Maps activity comment strings to active/inactive/unknown |
| `config/standard_text_manual_curation.csv` | 05 | Maps standard text values to active/inactive/unknown |

---

### Step 00 — Export ChEMBL activities (`00_export_chembl_activities.py`)

Exports the following raw tables from the local ChEMBL PostgreSQL database into `data/chembl_activities/`:

- `activities.csv`, `assays.csv`, `assay_parameters.csv`, `activity_stds_lookup.csv`
- `bioassay_ontology.csv`, `compound_structures.csv`, `docs.csv`
- `molecule_dictionary.csv`, `target_dictionary.csv`
- `target_components.csv`, `component_synonyms.csv`, `component_sequences.csv`, `source.csv`

For details on table contents see the [ChEMBL schema documentation](https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/schema_documentation.txt).

In addition, the script derives four frequency tables from the `ACTIVITIES` table and one reference table from the `ASSAYS` table, all saved to `data/chembl_processed/`:

1. `00_activity_comments.csv` — frequency of values in `activity_comment` (e.g. `"active"`)
2. `00_activity_std_units.csv` — frequency of (`standard_type`, `standard_units`) pairs (e.g. `"IC50 / nM"`)
3. `00_standard_units.csv` — frequency of `standard_units` values
4. `00_standard_text.csv` — frequency of values in `standard_text_value` (e.g. `"Compound metabolized"`)
5. `00_assay_descriptions.csv` — assay IDs mapped to their text descriptions and ChEMBL IDs (from `ASSAYS`)

The files `00_activity_comments.csv` and `00_standard_text.csv` were manually reviewed to assign an activity label to the most frequent entries, stored in `config/activity_comments_manual_curation.csv` and `config/standard_text_manual_curation.csv` (column `manual_curation_activity`): `1` = active, `-1` = inactive, `0` = inconclusive.

---

### Step 01 — Prepare manual curation files (`01_prepare_manual_files.py`)

Harmonizes activity type strings (uppercase, punctuation stripped) and maps units to their UCUM-compliant equivalents via `ucum_manual.csv`. Synonym activity types are collapsed using `synonyms.csv`. Produces two outputs in `data/chembl_processed/`:

- `01_activity_std_units_converted.csv` — (`activity_type`, `unit`) pairs with counts, ready for manual curation of biological direction
- `01_harmonized_types_map.csv` — each harmonized activity type mapped to the count and list of raw `standard_type` variants it collapses

> ⚠️ **Manual curation required before continuing.** The pipeline will exit with an error if the curated file is missing.

---

### Step 02 — Get compound info (`02_get_compound_info.py`)

Merges `compound_structures.csv` with `molecule_dictionary.csv` to map each `molregno` to its `chembl_id`. Calculates molecular weight from `canonical_smiles`. Output: `data/chembl_processed/02_compound_info.csv`.

---

### Step 03 — Standardize compounds (`03_standardize_compounds.py`)

Standardizes compound structures using the ChEMBL structure pipeline: canonicalizes SMILES, removes salts and solvents, and extracts the parent molecule. Molecular weight is recalculated from the standardized structure. Output: `data/chembl_processed/03_compound_info_standardized.csv`.

⏳ ETA: ~3–4 hours (main bottleneck of the pipeline).

---

### Step 04 — Merge activity and compounds (`04_merge_activity_and_compounds.py`)

Joins activity records with assay metadata, target information, and standardized compound data into a single unified table. Each row represents one experimental measurement with its associated compound, assay, and target fields. Output: `data/chembl_processed/04_activities_all_raw.csv`.

---

### Step 05 — Clean activities (`05_clean_activities.py`)

Produces the final preprocessed activity table `data/chembl_processed/05_activities_preprocessed.csv`. Performs the following steps:

1. **Filter invalid entries** — removes activities with missing canonical SMILES.
2. **Flag activity comments** — maps each `activity_comment` to active (1), inactive (-1), or unknown (0) using `config/activity_comments_manual_curation.csv`.
3. **Flag standard text** — same for `standard_text_value` using `config/standard_text_manual_curation.csv`.
4. **Convert units and values** — normalizes unit strings and converts raw values using formulas from `config/ucum_manual.csv`.
5. **Standardize relations** — maps `>=`, `>>`, `~` etc. to simplified `>`, `<`, `=`.
6. **Harmonize activity types** — strips punctuation/spaces and uppercases `standard_type`.
7. **Calculate pChEMBL** — recalculates pChEMBL values for records with `unit = umol.L-1` and a valid numeric value.
8. **Convert doc IDs** — replaces internal `doc_id` with `doc_chembl_id`.
9. **Map synonyms** — collapses activity type variants to their canonical name using `config/synonyms.csv`.
10. **Create text flag** — merges activity comment and standard text flags into a single `text_flag` column (1 / -1 / 0).

Also saves `data/chembl_processed/05_activity_std_units_curated_comments.csv`: a per-(`activity_type`, `unit`) summary of total counts and text-flagged records.

⏳ ETA: ~15 minutes.

---

### Step 06 — Calculate ECFPs (`06_calculate_ecfps.py`) *(optional)*

Computes ECFP fingerprints (radius 3, 2048 bits) for all standardized compounds using RDKit. Failed SMILES are skipped. Results are stored in HDF5 format. Output: `data/chembl_processed/06_chembl_ecfps.h5`.

⏳ ETA: ~15 minutes.

---

## Stage 2 — Generate pathogen datasets (run per pathogen)

```sh
conda activate camt
bash scripts/generate_datasets.sh --<pathogen_code>
# e.g. bash scripts/generate_datasets.sh --mtb
```

Supported codes are in `config/pathogens.csv`. All outputs go to `output/<pathogen_code>/`.

| Step | Script | What it does | Key output |
|------|--------|-------------|------------|
| 07 | `07_get_pathogen_assays.py` | Filters the full ChEMBL dataset to records matching the pathogen (by organism field and optional assay allowlist). ⏳ ~1 min | `07_chembl_raw_data.csv.gz`, `07_assays_raw.csv` |
| 08 | `08_clean_pathogen_activities.py` | Removes invalid records, assigns biological direction per (activity\_type, unit) pair, drops unmodelable activities. ⏳ ~5 min | `08_chembl_cleaned_data.csv.gz`, `08_assays_cleaned.csv` |
| 09 | `09_curate_assay_parameters.py` | **⚠️ GPU + ollama required.** Uses a local LLM to extract and standardize biological context for each unique assay: target type/name/ChEMBL ID, strain, ATCC ID, mutations, drug resistances, culture media. Processes assays in batches of 6. Already-processed assays are skipped on restart. Target ChEMBL IDs are resolved via a multi-strategy exact-match lookup (organism-prefix stripping, slash-split, gene-name reconstruction) — IDs are never guessed. | `09_assays_parameters_full.csv` |
| 10 | `10_calculate_assay_clusters.py` | Clusters compounds in each assay by ECFP4 similarity (BitBirch) at three Tanimoto thresholds (0.3, 0.6, 0.85). **Can run in parallel with step 09.** ⏳ ~10 min | `10_assays_clusters.csv` |
| 11 | `11_get_assay_overlap.py` | Computes pairwise compound overlap between assays with ≥ 50 compounds. ⏳ ~1 min | `11_assays_overlap.csv` |
| 12 | `12_prepare_assay_datasets.py` | Binarizes activities using expert cutoffs (`config/expert_cutoffs.csv`). Produces quantitative (`_qt`), qualitative (`_ql`), and mixed (`_mx`) datasets per assay. **⚠️ `expert_cutoffs.csv` must exist before running.** ⏳ ~5 min | `12_datasets.csv`, `datasets/datasets_qt.zip`, `datasets_ql.zip`, `datasets_mx.zip` |
| 13 | `13_lightmodel_individual.py` | Trains Random Forest classifiers (4-fold CV, Morgan fingerprints) on each binarized dataset. Evaluates AUROC. Adds ChEMBL decoys for active-enriched datasets. ⏳ ~30 min | `13_individual_LM.csv`, `13_reference_set.csv.gz`, `correlations/A/`, `correlations/B/` |
| 14 | `14_select_datasets_individual.py` | Selects the best cutoff per assay (AUROC > 0.7 threshold; prefers the mid cutoff). ⏳ < 1 min | `14_individual_selected_LM.csv` |
| 15 | `15_lightmodel_merged.py` | Merges assays that share experimental context (same activity type, unit, target, strain) and re-evaluates modelability for assays rejected in step 14. ⏳ ~10–30 min | `15_merged_LM.csv`, `correlations/M/`, `datasets/M/` |
| 16 | `16_select_datasets_merged.py` | Selects the best cutoff per merged group (same criteria as step 14). ⏳ < 1 min | `16_merged_selected_LM.csv` |
| 17 | `17_evaluate_correlations.py` | Computes pairwise similarity between all selected ORGANISM datasets using reference-set predictions, then greedily deduplicates redundant models. ⏳ ~5 min | `17_dataset_correlations.csv`, `17_final_datasets.csv` |
| 18 | `18_prepare_assay_master.py` | Assembles a master annotation table merging all per-step metadata, pipeline status flags, and selection traceability fields for every assay triplet. ⏳ < 1 min | `18_assays_master.csv` |
| 19 | `19_prepare_final_datasets.py` *(optional)* | Exports the selected datasets as simplified CSVs (SMILES + binary label only) in a single ZIP. ⏳ < 1 min | `19_final_datasets.zip` |
| 20 | `20_diagnosis.py` *(optional)* | Produces a diagnostic plot covering data quality, chemical diversity, assay coverage, and model correlations. ⏳ ~2 min | `20_diagnosis.png` |

### Config files used in Stage 2

| File | Step | Purpose |
|------|------|---------|
| `config/activity_std_units_manual_curation.csv` | 08 | Biological direction per (activity\_type, unit) pair |
| `config/expert_cutoffs.csv` | 12 | Binarization thresholds per (activity\_type, unit, target\_type, pathogen\_code) |
| `config/pathogens.csv` | all | Pathogen codes and names |

---

### Step 07 — Get pathogen assays (`07_get_pathogen_assays.py`)

Filters the full preprocessed ChEMBL dataset to extract all bioactivity records associated with the target pathogen. Matching is done by text search on `target_organism` and `assay_organism`, plus any assay ChEMBL IDs listed in `config/assays/<pathogen_code>.csv`.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `07_chembl_raw_data.csv.gz` | All activity records matching the pathogen. |
| `07_target_organism_counts.csv` | Frequency table of `target_organism` values. |
| `07_compound_counts.csv.gz` | Unique compounds with InChIKey, SMILES, and activity count. |
| `07_assays_raw.csv` | Per-(`assay_id`, `activity_type`, `unit`) summary with compound counts, text flags, and fraction of pathogen chemical space. |

⏳ ETA: ~1 minute.

---

### Step 08 — Clean pathogen activities (`08_clean_pathogen_activities.py`)

Applies sequential filters to the raw pathogen data and assigns biological direction per (`activity_type`, `unit`) pair. Activities are kept only if they have a valid direction or an active/inactive text flag.

Filters applied in order:

1. **Remove compounds with no SMILES** — activities without a valid structure are discarded.
2. **Remove empty activities** — activities lacking both a numeric value and a non-zero `text_flag` are discarded.
3. **Filter by consensus units** — only activities with units present in `data/chembl_processed/01_activity_std_units_converted.csv` (or no unit) are retained.
4. **Assign direction** — biological direction (-1 or +1) is assigned per (`activity_type`, `unit`) from `config/activity_std_units_manual_curation.csv`.
5. **Remove unmodelable activities** — activities with no direction and no active/inactive text flag are discarded.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `08_chembl_cleaned_data.csv.gz` | Cleaned activity records for the pathogen. |
| `08_activity_type_unit_comments.csv` | Per-(`activity_type`, `unit`) counts, cumulative proportion, text flag summary, and assigned direction. |
| `08_assays_cleaned.csv` | Same structure as `07_assays_raw.csv` with an added `direction` column, built on the cleaned data. |

⏳ ETA: ~5 minutes.

---

### Step 09 — Curate assay parameters (`09_curate_assay_parameters.py`)

> ⚠️ **Requires a GPU-enabled machine with [ollama](https://ollama.com/) running locally.**

Uses a local LLM to extract and standardize biological context for each (`assay_id`, `activity_type`, `unit`) triplet in `08_assays_cleaned.csv`. For each assay, a prompt is built from ChEMBL assay fields (`data/chembl_activities/assays.csv`) and publication metadata (`data/chembl_activities/docs.csv`), and the model is asked to return a strict JSON with the following fields:

| Field | Description |
|-------|-------------|
| `organism_curated` | Biological species explicitly stated (e.g. `Mycobacterium tuberculosis`) |
| `target_type_curated` | Curated target type: `SINGLE PROTEIN`, `ORGANISM`, `DISCARDED`, or the verbatim ChEMBL target type |
| `target_name_curated` | Target name explicitly stated in the annotations |
| `target_chembl_id_curated` | Target ChEMBL ID if explicitly stated |
| `strain` | Biological strain name (e.g. `H37Rv`); no catalog identifiers |
| `atcc_id` | ATCC identifier, formatted as `ATCC <number>` |
| `mutations` | Explicit mutations in one-letter format (e.g. `S450L`) |
| `known_drug_resistances` | Drugs for which resistance is explicitly stated |
| `media` | Growth or culture medium explicitly stated |

Target ChEMBL IDs extracted by the LLM are resolved against the ChEMBL target dictionary using a **multi-strategy exact-match lookup**: first by direct ID match, then by organism-prefix stripping, then by slash-split of compound names, and finally by gene-name reconstruction. IDs are never guessed — unresolvable names are left blank.

The script supports **resuming**: already-processed triplets in the output file are skipped. On LLM failure, an empty row is written and processing continues.

The LLM model is configured via `LLM_MODEL` in `src/default.py`.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `09_assays_parameters.csv` | Curated assay parameters for all (`assay_id`, `activity_type`, `unit`) triplets. |

⏳ ETA: ~5 seconds per assay on a GPU-enabled machine.

---

### Step 10 — Calculate assay clusters (`10_calculate_assay_clusters.py`)

For each (`assay_id`, `activity_type`, `unit`) triplet in `08_assays_cleaned.csv`, unique compound SMILES are retrieved from the cleaned activity data and clustered using ECFP4 fingerprints (2048 bits) with the BitBirch algorithm. The number of clusters is computed at three Tanimoto Coefficient thresholds (0.3, 0.6, 0.85), providing a measure of chemical diversity within each assay.

Steps 09 and 10 are **independent** and can be run in parallel after step 08.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `10_assays_clusters.csv` | Number of clusters per (`assay_id`, `activity_type`, `unit`) at each Tanimoto threshold. |

⏳ ETA: ~10 minutes.

---

### Step 11 — Get assay overlap (`11_get_assay_overlap.py`)

Computes pairwise compound overlap between assays in the cleaned pathogen dataset. Each assay is treated as an independent (`assay_id`, `activity_type`, `unit`) triplet. Only assays with 50 or more compounds are considered. For each pair, the number of shared compounds, a normalized overlap ratio, and whether both assays originate from the same publication are calculated.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `11_assays_overlap.csv` | Pairwise assay overlap table with compound counts, shared compounds, overlap ratio, and same-document flag. |

⏳ ETA: ~1 minute.

---

### Step 12 — Prepare assay datasets (`12_prepare_assay_datasets.py`)

Produces binarized, ML-ready compound-level datasets for each (`assay_id`, `activity_type`, `unit`) triplet. Requires the outputs of steps 08 and 09.

Each assay is processed through two parallel paths — **quantitative** (numeric values) and **qualitative** (text-based activity flags) — which are then combined into up to three dataset types per triplet.

#### Target type curation

The LLM-curated `target_type_curated` from step 09 is merged onto the cleaned assay table and then post-processed into a simplified `target_type_curated_extra` field, constrained to one of three values:

- `SINGLE PROTEIN` — also captures `PROTEIN COMPLEX` and `PROTEIN FAMILY`
- `ORGANISM` — whole-cell or phenotypic assays
- `DISCARDED` — assays that cannot be assigned a clear biological target

This simplified field is used as the key for expert cutoff lookup.

#### Expert cutoffs

Binarization thresholds are loaded from `config/expert_cutoffs.csv`, keyed by (`activity_type`, `unit`, `target_type_curated_extra`, `pathogen_code`). Multiple cutoffs per key are supported (semicolon-separated), generating one dataset per cutoff value. If no cutoff is defined for an assay, quantitative binarization is skipped.

> ⚠️ `expert_cutoffs.csv` must be manually created before running this step. Cutoffs are defined in the canonical unit space produced by step 05 (e.g. IC50 in `umol.L-1`). Because all values are normalized to a single unit per activity type in step 05, separate entries for nM, µM etc. are not needed.

#### Quantitative binarization

For assays with numeric values, a valid biological direction, and an expert cutoff:

1. **Adjust censored relations** — measurements on the wrong side of the direction (e.g. `IC50 > 100 µM` when lower = more active) are replaced with the assay extreme value and their relation set to `=`. This ensures they are correctly classified as inactive after binarization without distorting the compound selection step.
2. **Disambiguate compounds** — if a compound has multiple measurements, the most active one is kept: minimum for direction = −1, maximum for direction = +1.
3. **Binarize** — values are compared to the expert cutoff: `bin = 1` if `value ≤ cutoff` (direction = −1) or `value ≥ cutoff` (direction = +1).

#### Qualitative binarization

Built from the `text_flag` column (derived from activity comments and standard text values in step 05). Aggregated to compound level:
- If a compound has conflicting labels (both active and inactive), a `ValueError` is raised.
- Label priority: `1 > −1 > 0`. Compounds with only neutral flags (0) are removed.

#### Dataset types

| Type | Condition | File suffix |
|------|-----------|-------------|
| `quantitative` | Numeric values + valid direction + expert cutoff, no qualitative data | `_qt_<cutoff>.csv.gz` |
| `qualitative` | Text flags only, or no cutoff/direction available | `_ql.csv.gz` |
| `mixed` | Both quantitative and qualitative data available | `_mx_<cutoff>.csv.gz` |
| `none` | Quantitative data exists but no cutoff or direction; no qualitative data | — |

For **mixed** datasets, qualitative inactives not already present in the quantitative set are appended to it. Only inactives are added this way — qualitative actives without a numeric measurement are excluded to avoid inflating the active class with unquantified data.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `12_datasets.csv` | One row per (assay, cutoff) — dataset type, compound counts, positive counts, ratios. |
| `12_assay_data_info.csv` | One row per assay — dataset type, relation counts (`=`, `<`, `>`), value distribution percentiles. |
| `datasets/datasets_qt.zip` | All quantitative dataset files. |
| `datasets/datasets_ql.zip` | All qualitative dataset files. |
| `datasets/datasets_mx.zip` | All mixed dataset files. |

⏳ ETA: ~5 minutes.

---

### Step 13 — Individual light modeling (`13_lightmodel_individual.py`)

Evaluates the modelability of each binarized assay dataset from step 12 using Random Forest classifiers trained on Morgan (ECFP) fingerprints. Produces per-assay AUROC scores and reference set predictions used for downstream correlation analysis.

#### Modeling conditions

Datasets are evaluated under two independent conditions that differ in size requirements and class balance strategy:

| | Condition A | Condition B |
|--|-------------|-------------|
| Dataset types | quantitative or mixed | quantitative or mixed |
| Min compounds | ≥ 1,000 | any |
| Min actives | ≥ 50 | ≥ 100 |
| Active ratio | 0.001–0.5 | ≥ 0.5 |
| Decoys added | No | Yes |
| Purpose | Large, balanced datasets | Active-enriched datasets, decoy-balanced |

The same assay can qualify under both conditions and will be modeled independently in each.

#### Condition B decoy strategy

For datasets where actives are over-represented (ratio ≥ 0.5), random ChEMBL compounds that have **never been tested against the pathogen** are added as decoys until the active ratio reaches 10% (`decoy_ratio = 0.1`). This prevents the model from learning a biased baseline while keeping the active set intact.

#### Per-assay workflow

For each qualifying dataset:
1. Load the quantitative or mixed dataset from the step 12 zip archives
2. *(Condition B only)* Sample and append decoys from the pathogen-naive ChEMBL pool
3. Run **4-fold stratified cross-validation** with Random Forest (100 trees) → mean AUROC ± std
4. Train a **final model on all data** → predict activity probabilities on the reference set

#### Reference set

The reference set is a fixed sample of **10,000 compounds that have never been screened against the pathogen**, drawn from all ChEMBL compounds with available fingerprints. Using pathogen-naive compounds avoids data leakage: correlations between model predictions reflect shared biological signal rather than shared training data. The set is sampled with a fixed seed (42) for reproducibility.

> The alternative of using pathogen-tested compounds was considered but rejected: many of those compounds appear in the training sets of the models being compared, which would make correlations partly reflect memorization rather than generalization.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `13_individual_LM.csv` | One row per (assay, condition) — AUROC mean and std, dataset metadata. |
| `13_reference_set.csv.gz` | The 10,000 pathogen-naive reference compounds. |
| `correlations/A/` | Reference set prediction probabilities (`.npz`) for condition A models. |
| `correlations/B/` | Reference set prediction probabilities (`.npz`) for condition B models. |

⏳ ETA: ~30 minutes.

---

### Step 14 — Select individual datasets (`14_select_datasets_individual.py`)

Selects the single best dataset per assay triplet from the step 13 results, applying an AUROC threshold and a cutoff preference rule.

#### Selection criteria

Only assays with a best AUROC > 0.7 (across all cutoffs tested) are retained. For assays with multiple expert cutoffs, one cutoff is chosen per the following rule:

- **Prefer the mid cutoff** (the second value in the `expert_cutoffs.csv` list for that assay) as a default, since it represents the central activity threshold and is more interpretable.
- **Switch to the best cutoff** only if the mid cutoff is unavailable or if the best cutoff achieves an AUROC more than 0.1 higher than the mid cutoff — indicating a meaningfully better model at that threshold.

The flag `is_mid_cutoff` records which rule was applied for each selected dataset.

#### Per-label independence

Selection is performed independently for conditions A and B. The same assay triplet can therefore appear in both — with potentially different cutoffs and AUROC scores — since the two conditions model different compound sets (condition B adds decoys).

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `14_individual_selected_LM.csv` | Selected datasets with cutoff, AUROC, and dataset statistics. |

⏳ ETA: < 1 minute.

---

### Step 15 — Merged light modeling (`15_lightmodel_merged.py`)

Attempts to rescue assay triplets that were not accepted in step 14 by merging them into larger combined datasets and re-evaluating their modelability.

Many assays in ChEMBL are too small individually to meet the compound and positive thresholds of conditions A or B. However, assays that share the same experimental context (activity type, unit, direction, assay type, target type, BAO label, and strain) can be combined into a single dataset. If the union of their compounds is large enough, a viable classifier can still be built.

#### Merge candidate identification

Assay triplets that were **not** accepted in step 14 are grouped by shared experimental metadata. Two grouping strategies are applied independently:

- **ORGANISM** — groups by (`activity_type`, `unit`, `direction`, `assay_type`, `target_type_curated_extra`, `bao_label`, `strain`). Captures whole-cell / phenotypic assay merges.
- **SINGLE PROTEIN** — same keys plus `target_chembl_id`. Ensures that only assays against the same protein target are combined.

A group qualifies for merging if the **union** of compounds across all member assays exceeds 1,000 and the group contains at least 2 assays. Groups are ranked by union size.

#### Merging and binarization

For each qualifying group and each applicable expert cutoff:

1. Quantitative and mixed datasets for all member assays are loaded from the step 12 zip archives and concatenated.
2. Across the merged quantitative records, **the most active measurement is kept per compound** (min for direction = −1, max for direction = +1).
3. Qualitative inactives from mixed datasets are appended, with compound-level deduplication (existing records take priority).
4. Datasets with fewer than 50 positives after merging are skipped.

#### Decoy strategy

If the active ratio in the merged dataset exceeds 50%, random ChEMBL decoys (pathogen-naive compounds) are added until the ratio reaches 10% — the same strategy as condition B in step 13.

#### Per-group workflow

For each qualifying group:
1. Build the merged dataset as described above
2. *(If active ratio > 0.5)* Sample and append decoys
3. Run **4-fold stratified cross-validation** with Random Forest → mean AUROC ± std
4. Train a **final model on all data** → predict activity probabilities on the reference set (loaded from `13_reference_set.csv.gz`)
5. Save merged dataset and reference predictions

Merged datasets are named `M_ORG<i>` (ORGANISM) or `M_SP<i>` (SINGLE PROTEIN), suffixed with the expert cutoff value.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `15_merged_LM.csv` | One row per merged group × cutoff — AUROC, compound counts, assay keys, metadata. |
| `correlations/M/` | Reference set prediction probabilities (`.npz`) for each merged model. |
| `datasets/M/` | Merged dataset files (`.csv.gz`), one per group × cutoff. |

⏳ ETA: ~10–30 minutes depending on the number of qualifying groups.

---

### Step 16 — Select merged datasets (`16_select_datasets_merged.py`)

Selects the single best dataset per merged group from the step 15 results, applying the same AUROC threshold and mid-cutoff preference rule used in step 14 for individual datasets.

#### Selection criteria

Only merged groups with a best AUROC > 0.7 (across all cutoffs tested) are retained. For groups with multiple expert cutoffs, one cutoff is chosen per the following rule:

- **Prefer the mid cutoff** (the second value in the `expert_cutoffs.csv` list) as a default, since it represents the central activity threshold and is more interpretable.
- **Switch to the best cutoff** only if the mid cutoff is unavailable or if the best cutoff achieves an AUROC more than 0.1 higher — indicating a meaningfully better model at that threshold.

The flag `is_mid_cutoff` records which rule was applied for each selected dataset.

#### Coverage tracking

The script builds an assay → compound mapping from `08_chembl_cleaned_data.csv.gz` and reports the chemical space coverage before and after selection, separately for ORGANISM and SINGLE PROTEIN groups.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `16_merged_selected_LM.csv` | Selected merged datasets with cutoff, AUROC, and group metadata. |

⏳ ETA: < 1 minute.

---

### Step 17 — Evaluate dataset correlations (`17_evaluate_correlations.py`)

Evaluates pairwise similarity between all selected ORGANISM datasets (individual A/B from step 14 and merged M from step 16), then applies a greedy deduplication to select a non-redundant final set.

#### Why only ORGANISM assays?

Phenotypic whole-organism assays are the primary endpoint of interest for antimicrobial drug discovery. SINGLE PROTEIN datasets are modeled in earlier steps but excluded from the final selection here. They remain accessible via the master table (step 18) with their pipeline status annotated.

#### Similarity metrics

Similarity between any two models is measured using their predictions on the **reference set** — 10,000 pathogen-naive compounds saved in step 13. Three complementary metrics are computed:

| Metric | Description |
|--------|-------------|
| `spearman` | Rank correlation of predicted probabilities across all 10,000 reference compounds. |
| `hit_overlap_1000` | Normalised overlap of the top-1,000 predictions above chance (broad hit agreement). |
| `hit_overlap_100` | Normalised overlap of the top-100 predictions above chance (high-confidence hit agreement). |
| `compound_overlap` | Fraction of training compounds shared between the two datasets (Jaccard-like, normalised by the smaller set). |

All pairwise combinations are computed and saved to `17_dataset_correlations.csv`.

#### Greedy deduplication

Datasets are sorted by size (largest first) and processed greedily: a dataset is discarded if it is simultaneously **model-correlated** and **compound-overlapping** with an already-selected dataset:

```
(spearman + hit_overlap_1000 + hit_overlap_100) / 3 > 0.5
AND compound_overlap > 0.5
```

Using `AND` means that high model correlation alone is not enough to discard a dataset — there must also be substantial compound overlap. Two models that are highly correlated but trained on different compounds may be capturing the same biology across different chemical series, which is informative rather than redundant.

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `17_dataset_correlations.csv` | All pairwise similarity scores between selected ORGANISM datasets. |
| `17_final_datasets.csv` | All ORGANISM datasets with a `selected` flag indicating which passed deduplication. |

⏳ ETA: ~5 minutes (scales quadratically with the number of datasets).

---

### Step 18 — Prepare assay master table (`18_prepare_assay_master.py`)

Assembles a comprehensive per-assay annotation table that integrates all metadata and pipeline status flags across every step. This is the primary reference output for understanding what happened to each assay in the pipeline.

#### Contents

The master table contains one row per (`assay_id`, `activity_type`, `unit`) triplet and merges the following sources:

| Source file | Fields added |
|-------------|--------------|
| `08_assays_cleaned.csv` | Base assay metadata: assay type, target, direction, compound counts, text flags |
| `09_assays_parameters.csv` | LLM-curated fields: organism, target type, strain, ATCC ID, mutations, media |
| `10_assays_clusters.csv` | Chemical diversity: cluster counts at Tanimoto 0.3, 0.6, 0.85 |
| `12_assay_data_info.csv` | Dataset type, relation counts, activity value distribution (percentiles) |

An `evaluated_cutoffs` column is added from `config/expert_cutoffs.csv`, listing the cutoffs that were tested for each assay.

#### Pipeline status annotations

Three `comment_A`, `comment_B`, `comment_M` columns annotate each assay's journey through the pipeline under conditions A, B, and M respectively. Each takes one of:

| Status | Meaning |
|--------|---------|
| `Not considered` | Did not meet the size/balance criteria for modeling |
| `Considered but not selected` | Modeled but AUROC ≤ 0.7 |
| `Considered and selected, but not in correlation analysis (non-ORGANISM)` | AUROC > 0.7 but excluded from step 17 (SINGLE PROTEIN) |
| `Considered and selected, but discarded due to high correlation` | Good model but redundant with a larger dataset |
| `Considered, selected, and retained in final selection` | Part of the final non-redundant dataset set |

Outputs are saved to `output/<pathogen_code>/`:

| File | Description |
|------|-------------|
| `18_assays_master.csv` | One row per assay triplet with all metadata and pipeline status columns. |

⏳ ETA: < 1 minute.

---

### Step 19 — Prepare final datasets (`19_prepare_final_datasets.py`) *(optional)*

Exports all selected datasets as simplified CSVs containing only SMILES and binary activity labels, bundled into a single ZIP file. Intended for researchers who want ML-ready data without navigating the full pipeline outputs.

Output: `output/<pathogen_code>/19_final_datasets.zip`

⏳ ETA: < 1 minute.

---

### Step 20 — Diagnosis (`20_diagnosis.py`) *(optional)*

Produces a diagnostic plot summarizing data quality, chemical diversity, assay coverage, and model correlations for the pathogen dataset.

Output: `output/<pathogen_code>/20_diagnosis.png`

⏳ ETA: ~2 minutes.

---

## Accessing ML-ready datasets

After the full pipeline, researchers can access datasets without manually sifting files:

```sh
# Recommended: simplified export (SMILES + binary label)
python scripts/19_prepare_final_datasets.py <pathogen_code>
unzip output/<pathogen_code>/19_final_datasets.zip
```

Finally, script 21 is required to generate files for the pubchem pipeline, but all pathogens need to be processed minimum to step 07 for it to work.

For programmatic access with full metadata, see `docs/dataset_access_guide.md`.
