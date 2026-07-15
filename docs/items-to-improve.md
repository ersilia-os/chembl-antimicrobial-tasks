# Items to improve

A running list of data-quality / pipeline issues found during audit. Each entry
records what was observed, how it was measured, its impact, and open decisions.
Add new items as they are found; do not remove entries without sign-off.

---

## 1. LLM over-assigns `assay_organism_curated = "Homo sapiens"`

**Found:** 2026-07-15 (while auditing `28_trace_sample_datasets.py`).

**What it is.** Step 09 (`09_curate_assay_parameters.py`) asks the LLM to relabel
the assay organism as `Homo sapiens` when a description points to a human-cell
cytotoxicity / counter-screen. This instruction over-fires: it relabels many
assays that are plainly antibacterial/antifungal/antiparasitic tests *against the
pathogen itself*.

**Measurement.** Across all 16 pathogens (`21_curation_summary.csv`), 593 datasets
have `assay_organism` = the pathogen but `assay_organism_curated` = `Homo sapiens`.
Splitting by whether the description contains any human/cytotoxicity cue
(`human|hela|hepg2|vero|fibroblast|cytotox|mammalian|cell line|macrophage|…`):

| group | count |
|---|---|
| total pathogen → Homo sapiens | 593 |
| description has a human/cytotox cue (plausibly correct) | 121 |
| no human cue (suspicious) | 472 |
| no human cue **and** description names the pathogen (almost certainly wrong) | 458 |

By pathogen: ecoli 169, saureus 126, mtuberculosis 85, pfalciparum 78,
paeruginosa 41, calbicans 39, kpneumoniae 21, abaumannii 19, others ≤6.

**Representative errors** (all are pathogen assays relabelled `Homo sapiens`):

| assay | description (abridged) | fate |
|---|---|---|
| `CHEMBL2044012` | Antibacterial activity against *Klebsiella pneumoniae* … agar well diffusion | retained (ORGANISM) |
| `CHEMBL3875364` | Antibacterial activity against *Streptococcus pneumoniae* RCMB 010010 … agar well diffusion | retained (ORGANISM) |
| `CHEMBL3767960` | Antibacterial activity against *Staphylococcus aureus* ATCC 25923 … microdilution | retained (ORGANISM) |
| `CHEMBL3766769` | Antitubercular activity against *Mycobacterium bovis* BCG ATCC 35743 … NR assay | retained (ORGANISM) |
| `CHEMBL1110336` | Antifungal activity against *Candida albicans* by XTT assay | discarded (≤5 mol) |
| `CHEMBL871974` | Parasite growth inhibition against *Plasmodium falciparum* GP1 | discarded (≤5 mol) |
| `CHEMBL5373389` | Resistance development in *Pseudomonas aeruginosa* PAO1 … FoR assay | discarded (qual-only) |
| `CHEMBL1047463` | Binding affinity to *E. coli* pUC19 DNA … agarose gel | discarded (not ORGANISM) |

**Impact — real but currently inert.** `assay_organism_curated` is never used as a
classifier or filter anywhere in the pipeline (verified: it appears only in
column-carry lists in steps 16/19/21, never in a groupby / filter / comparison /
merge key). Inclusion/exclusion is decided by `target_type` → `target_type_curated`
→ `target_type_curated_extra` (must be `ORGANISM`); strain merging uses
`assay_strain_norm`. So a mislabelled assay still enters the correct pathogen
dataset. The mislabel corrupts a **display/provenance** column (audit panel, assay
master), not the built datasets — it is a credibility issue, not data leakage.

**Open decisions (not actioned):**
- Error definition / scope: is "raw names the pathogen + no human cue" the right
  rule? Should the ~458 be dumped to a reviewable CSV?
- Remedy: (a) leave as-is (inert); (b) flag likely-mislabels in the audit panel;
  (c) tighten the step-09 human-cell prompt — this would require re-running step 09.

---

## 2. `redundant` leftovers are dropped without reconciling their labels

**Found:** 2026-07-15 (auditing step-25 leftover handling; see
[`2026-07-15_redundant-conflict-diagnostics.md`](2026-07-15_redundant-conflict-diagnostics.md),
tool [`tools/redundant_conflict_audit.py`](../tools/redundant_conflict_audit.py)).

**What it is.** Step 25 drops a leftover as `redundant` when all its molecules are already in a
modelled pool (`n_residual == 0`), keeping only the pool's existing label. That assay's own reading
is never folded in, so any disagreement is discarded.

**Measurement.** Over 2,455 redundant leftovers / 29,811 matched molecules: raw consensus conflict
**38.4 %** (strict 27.7 %, any 54.4 %). Pipeline-faithful (merge + re-binarise, every containing
pool, 66,902 molecule × pool): **flip_vs_delivered 17.3 %**, of which the dropped assay's pure
marginal effect is **9.7 %**; flips are overwhelmingly **→ inactive** (10,408 vs 1,176). Driver is
numeric **MIC-vs-MIC** disagreement at the 10 µM DR cutoff (MIC = 21k of 30k matched, 41 %), not
qualitative assays.

**Impact.** ~10 % of these "already covered" molecules would take a different (usually inactive)
label if the redundant assay were merged in. Whether that is signal to recover or a small assay
correctly overruled by the aggregate is a scientific judgement.

**Open decisions (not actioned):**
- (a) leave as-is; (b) fold redundant-leftover labels back via `merge_datasets` majority vote and
  measure CV-AUROC impact on the dev organisms before adopting; (c) treat the → inactive skew as a
  cutoff-sensitivity signal and revisit the DR cutoff instead.

## 3. Step-25 "base wins / never re-average" pool-growth policy

**Found:** 2026-07-15 (same audit).

**What it is.** When step 25 merges a leftover into a pool it appends only the leftover's *residual*
molecules; for a molecule already in the pool it keeps the base value and discards the leftover's
reading (`drop_duplicates("inchikey")`), even though the leftover is listed as a full pool member in
`25_pool_members.csv`. A clean re-merge of the listed members can therefore differ from the shipped
label.

**Measurement.** Pool-wide, de-duplicated: only **267 / 697,471 molecule × pool = 0.04 %** of
shipped labels would flip if step 25 re-averaged (232 distinct molecules, ≤ 0.3 % per pathogen). The
larger 8.1 % seen per-leftover in the merge-flip analysis is multiplicity — the affected molecules
are promiscuous (~42 redundant leftovers each).

**Impact — negligible.** The policy is effectively inert at the pool level (0.04 %). Recorded for
completeness / future editions.

**Open decisions (not actioned):**
- (a) leave as-is (recommended given 0.04 %); (b) re-average existing molecules for cleanliness if
  step 25 is ever refactored — not worth a dedicated re-run on its own.

## 4. Reconsider the `MIN_COMPOUNDS = 100` floor in the PubChem pipeline

**Found:** 2026-07-15 (tracing calbicans `CHEMBL1794297` in the dataset-trace panel).

**What it is.** The companion PubChem pipeline
([`pubchem-antimicrobial-tasks`](https://github.com/ersilia-os/pubchem-antimicrobial-tasks),
`02_bioassays_not_in_chembl.py`) defers a ChEMBL-linked assay to PubChem only when PubChem has
**both** >10 % more compounds (`CHEMBL_MISMATCH_THRESHOLD = 0.1`) **and** at least
`MIN_COMPOUNDS = 100` compounds. An assay whose PubChem counterpart is larger than ChEMBL's but
still below 100 compounds clears the mismatch test yet is dropped by the size floor — so it is never
written to `config/pubchem_aids/chembl_assays_in_pubchem_<code>.csv`, and this repo keeps its smaller
ChEMBL slice instead of the fuller PubChem data.

**Example.** calbicans `CHEMBL1794297` (AC50, DR) is a ChEMBL import (`src_assay_id`) of PubChem
**AID 488930** — the Broad Institute "Reversing Antifungal Drug Resistance" fluconazole
re-sensitisation screen. ChEMBL holds **13** compounds; PubChem holds **34** (13 active). The
mismatch is flagged (161 % surplus, in `chembl_cpds_mismatch_candida_albicans.csv`) but then dropped
because 34 < 100, so ChEMBL's 13-molecule version is what ends up in the pool.

**Impact.** Sub-100 PubChem assays that materially outsize their ChEMBL counterpart are silently kept
at the smaller ChEMBL size rather than replaced by the fuller PubChem data. How many assays this
touches, and how many extra molecules it would recover, has not been quantified.

**Open decisions (not actioned):**
- Lowering `MIN_COMPOUNDS` (or applying a lower floor specifically to the compound-mismatch branch)
  would let these larger PubChem versions be delivered. It is a **threshold decision affecting both
  repos** and shifts many assays across the ChEMBL/PubChem boundary — quantify the effect on both
  pipelines and get sign-off before changing it.
