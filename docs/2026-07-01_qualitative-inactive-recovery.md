# Recovering qualitative-inactive molecules (steps 22–23)

**Date:** 2026-07-01
**Scope:** `scripts/22_build_curated_datasets.py`, `scripts/23_binarise_datasets.py`,
`src/dataset_utils.py` (`binarize_conservative`), `src/plot_utils.py` (`curation_stats`).

## Problem

Steps 22–23 build one dataset per retained `(assay_chembl_id, activity_type)` pair, keeping
**only numeric records** (`num = acts[acts["is_numeric"]]`) and deduplicating by InChIKey. A
molecule whose records in that pair were **entirely qualitative** (a `text_flag` active/inactive
call with no measured value) was therefore dropped silently — it never entered the written
dataset, and nothing in the curation summary recorded the loss.

Audit of the drop (before this change, all 15 pathogens):

- **3,118** text-only molecule-slots dropped from retained pairs — 0.25 % of ~1.22 M slots.
- The loss was **directional**: **98.3 % carried a qualitative *inactive* flag** (`text_flag = -1`),
  only 52 were qualitative actives.
- It landed on **small datasets** (median 13 molecules), so per-dataset it was not noise.
- Net effect: because binarisation happens *after* the drop, removing negatives **inflated the
  active ratio** and turned a number of genuinely mixed datasets into apparent all-active /
  degenerate ones (e.g. `saureus/CHEMBL3421263/GI`: 10/10 active, but 21 qualitative inactives
  had been dropped).

## Decision (curator sign-off)

For each **retained** pair, a molecule (InChIKey) with **no numeric record** is:

- **kept as INACTIVE** iff its qualitative records are purely inactive (`text_flag == -1`,
  no active flag);
- **discarded** if it carries **any** active flag (`text_flag == 1`) — it cannot be placed on
  the active side without a measured value (this includes the rare active+inactive conflict);
- left dropped if it has no ±1 flag (no signal).

Recovered inactives **count toward** the `MIN_MOLECULES = 5` size threshold (filter 6), so a pair
with few numeric records but enough qualitative inactives can now qualify.

## Implementation

- **Step 22** computes, per `(pair, InChIKey)`, `has_numeric / has_active_flag / has_inactive_flag`;
  derives the recovered set; changes filter 6 to test `n_mol_dataset = n_mol_numeric +
  n_recovered_inactive`; and appends recovered rows (`value`/`relation = NaN`, pair `unit` /
  `direction` / `category`) to each written dataset. New summary columns:
  `n_recovered_inactive`, `n_text_only_active_discarded`, `n_mol_dataset`.
- **Step 23** binarises unchanged; `binarize_conservative` now labels value-less rows inactive
  explicitly (`& value.notna()`). New summary column: `n_qual_inactive`.
- **`curation_stats`** reports retained sizes from `n_mol_dataset`.

## Effect after re-running all 15 pathogens

Verification invariants (32,557 binarised files): **0** value-less rows labelled active; **0**
`n_molecules` vs `n_mol_dataset` mismatches.

| quantity | before | after | Δ |
|---|---|---|---|
| binarised DR/SP datasets | 32,311 | 32,557 | +246 |
| degenerate datasets | 11,325 | 11,383 | +58 |
| total molecule-slots (DR/SP) | 1,219,976 | 1,225,920 | +5,944 |
| active labels | 282,740 | 283,000 | +260 |

Qualitative inactives recovered: **5,145** across retained pairs (**5,130** land in binarised
DR/SP datasets). Qualitative-active text-only molecules discarded: **24,512**.

Decomposition of the degenerate change (14 pathogens with exact before/after snapshots):

- **51 existing datasets flipped degenerate → non-degenerate** (recovered inactives supplied a
  negative class) — the intended de-inflation.
- **0 datasets flipped non-degenerate → degenerate** (adding inactives can never create an
  all-active set).
- **246 newly-retained datasets** (cleared ≥5 only once inactives count): **135 usable**
  (both classes) + 111 genuinely all-inactive screens.

So the headline "+58 degenerate" is not a regression: it is 111 previously-discarded,
qualitative-inactive-heavy screens now large enough to keep, partly offset by 51 datasets that
were un-inflated. The active label count is essentially unchanged for pre-existing datasets
(the +260 comes entirely from newly-retained pairs that contain measured actives).

## Note on the earlier audit figure

The pre-change per-dataset audit table labelled *all* text-only drops "dropped_inactives"; the
strict rule here recovers only the purely-inactive subset. `saureus/CHEMBL3421263/GI` and
`smansoni/CHEMBL5331366/ACTIVITY`, cited earlier as candidates to un-degenerate, turned out to be
qualitative **actives** and are correctly discarded — they remain degenerate.
