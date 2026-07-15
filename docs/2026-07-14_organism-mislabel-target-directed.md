# ChEMBL ORGANISM assays that are actually target/pathway-directed

**Status:** **resolved for the 6 confirmed cases** (2026-07-14). A config-driven per-assay target-type
override was added; the 5 target-directed assays are relabelled SINGLE PROTEIN and the 1 whole-cell
counterscreen is excluded. The broader detection gap (below) remains open — no automated sweep was run.

## The example — CHEMBL1794497 (E. coli)

- **Description:** `PUBCHEM_BIOASSAY: High Throughput Screen for Tat Transport Inhibitors Measured in
  Microorganism System Using Plate Reader ...` (related PubChem AID488938).
- **ChEMBL structured fields:** `assay_type=F`, `target_type=ORGANISM`,
  `target_chembl_id=CHEMBL354` (E. coli organism), `bao_label=organism-based format`,
  `activity_type=AC50`, `unit=umol.L-1`, 189 compounds.
- **Problem:** it is a **pathway-directed** reporter screen (Tat transport inhibition), but ChEMBL
  assigns the target to the whole organism, so our pipeline ingests it as generic antibacterial
  activity and pools it with phenotypic MIC-type data.

## Why our curation does not catch it

Two independent gates both pass it through unchanged:

1. `scripts/09_curate_assay_parameters.py` → `classify_assay`: any assay with `target_type=="ORGANISM"`
   and a matching organism is routed to the `text_only` path, where `prefill_target_fields` hard-sets
   `target_type_curated="ORGANISM"`. **The LLM never re-reads the description to reconsider the target
   type for ORGANISM assays.**
2. `src/pathogen_utils.py` → `extra_curation_target_type`: the `ORGANISM` branch returns `'ORGANISM'`
   unconditionally, so even an LLM correction could not override ChEMBL's raw label.

## How common is it — and can we detect it automatically?

Scanned all **93,014 unique ORGANISM assays** across the 15 pathogens. There is **no reliable
structured flag**: CHEMBL1794497's structured fields are all "normal" for an organism assay; its
target-directed nature lives only in the free-text description.

| Candidate signal | Count | Verdict |
|---|---|---|
| `assay_type=B` (Binding) + ORGANISM | 211 | Noisy — mostly genuine whole-cell assays (MIC, membrane-permeabilization) tagged "Binding". |
| BAO `single protein format` + ORGANISM | 1 | A plain P. falciparum IC50; not target-directed. |
| Molecular affinity types (Ki/Kd/Kb) + ORGANISM | 0 | ChEMBL routes these to SINGLE PROTEIN, so none survive as ORGANISM. |
| `activity_type=AC50` + ORGANISM (hypothesis: AC50 ⇒ enzyme) | 14 | **Refuted.** 13/14 are whole-cell (antifungal-resistance/toxicity screens, growth counterscreen, anti-adhesion, time-kill); only CHEMBL1794497 is target-directed. AC50 in PubChem HTS deposits is a generic dose-response label, not enzyme-specific. ~7% precision. |
| Text `"...Measured in Microorganism/Cell System"` (this PubChem family) | 12 | Heterogeneous — mixes real phenotypic assays with counterscreens and the Tat one. |

**Units are not a usable discriminator:** the affinity types that would be diagnostic (Ki/Kd/Kb) never
reach ORGANISM; the concentration types that do (IC50, EC50, ED50, AC50) are all used for whole-cell
readouts.

**Conclusion:** comprehensive detection would require a text/LLM pass over ORGANISM descriptions —
exactly the step the pipeline currently shortcuts. Keyword sweeps flag thousands with heavy false
positives.

## Resolution applied (2026-07-14, the 6 confirmed cases)

**Relabelled ORGANISM → SINGLE PROTEIN** (target-directed; dropped from Stage-4 ORGANISM datasets, recorded
as discard-step 1 "not ORGANISM (SINGLE PROTEIN)"): `CHEMBL1638137`, `CHEMBL1673903`, `CHEMBL1673904`
(pfalciparum), `CHEMBL2434716` (kpneumoniae), `CHEMBL1794497` (ecoli, Tat).

**Excluded entirely** (whole-cell counterscreen): `CHEMBL1794494` (ecoli).

Mechanism:
- New `config/target_type_overrides.csv` (`assay_id,target_type`) → loaded by
  `load_target_type_overrides()` and consulted at the top of `extra_curation_target_type()`
  (`src/pathogen_utils.py`), wired into both live recompute sites (`13_prepare_assay_datasets.py`,
  `21_build_curated_datasets.py`). The ORGANISM hard-lock is now overridable per assay.
- Counterscreen via `config/assays_to_exclude/ecoli.csv` (step 07). Its row was also removed from
  ecoli's `09_assays_parameters_full.csv` so Stage 3's bijective join (`add_target_type_curated`) stays
  consistent without a step-09 (GPU) re-run.
- Regenerated Stage 4 (21–28) for ecoli/pfalciparum/kpneumoniae; no step-09 re-run needed.

## Open decisions (still deferred)

- Detection scope for the *rest* of the corpus: LLM re-read of ORGANISM descriptions vs. keyword candidate
  list vs. per-assay manual. Only the 6 cases above are handled; no exhaustive sweep was run.
