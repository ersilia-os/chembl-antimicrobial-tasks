# Validation report - pipeline steps 21 & 22

Generated: 2026-07-09  
Tool: `tools/validate_datasets.py`  
Pathogens: abaumannii, calbicans, campylobacter, ecoli, efaecium, enterobacter, hpylori, kpneumoniae, mtuberculosis, ngonorrhoeae, paeruginosa, pfalciparum, saureus, smansoni, spneumoniae

**Overall: ALL CHECKS PASS**

## abaumannii

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 2991 pairs = 855 retained + 2136 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## calbicans

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 11051 pairs = 2835 retained + 8216 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## campylobacter

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 430 pairs = 47 retained + 383 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 47) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## ecoli

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 21484 pairs = 5604 retained + 15880 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## efaecium

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 1880 pairs = 650 retained + 1230 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## enterobacter

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 1724 pairs = 473 retained + 1251 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## hpylori

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 693 pairs = 129 retained + 564 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## kpneumoniae

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 5404 pairs = 1862 retained + 3542 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## mtuberculosis

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 10436 pairs = 3258 retained + 7178 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## ngonorrhoeae

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 511 pairs = 50 retained + 461 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 50) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## paeruginosa

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 12362 pairs = 3529 retained + 8833 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## pfalciparum

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 9535 pairs = 2921 retained + 6614 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## saureus

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 29871 pairs = 8783 retained + 21088 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## smansoni

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 809 pairs = 92 retained + 717 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## spneumoniae

| Check | Title | Result | Violations | Note |
|---|---|---|---|---|
| A1 | accounting (retained XOR reason; counts sum) | PASS | 0 | 4150 pairs = 1469 retained + 2681 discarded |
| A2 | coverage summary<->files | PASS | 0 | missing files: 0, orphan files: 0 |
| A3 | attrition monotonic non-increasing | PASS | 0 |  |
| B1 | InChIKey unique per dataset | PASS | 0 |  |
| B2 | row count == n_mol_dataset | PASS | 0 |  |
| B3 | no cross-assay contamination | PASS | 0 |  |
| B5 | no nulls in inchikey/smiles (+direction for DR/SP) | PASS | 0 |  |
| B6 | qualifiers valid ({=,>,<} or NaN) | PASS | 0 |  |
| B4 | retained DR/SP pairs have scalar unit & direction | PASS | 0 |  |
| D | no retained pair violates a curation filter | PASS | 0 |  |
| C | dedup keeps most-active (sample of 60) | PASS | 0 |  |
| E3 | coverage retained DR/SP <-> binarised files | PASS | 0 | missing: 0, orphan: 0 |
| E1 | binarisation drops no rows | PASS | 0 |  |
| E2 | molecule<->value identity preserved | PASS | 0 |  |
| F1 | bin recomputed (conservative) matches stored | PASS | 0 |  |
| F5 | no unproven actives (prove-the-hit) | PASS | 0 |  |
| F4 | bin in {0,1} | PASS | 0 |  |
| F2 | summary cutoff == config cutoff | PASS | 0 |  |
| F3 | summary counts/ratio/degenerate consistent | PASS | 0 |  |
| G | conservative labelling (DR): actives all proven | PASS | 0 |  |
| G | conservative labelling (SP): actives all proven | PASS | 0 |  |

## Worked examples

### abaumannii - step 21 dedup keeps most-active record

```
pair: CHEMBL1640427 / MIC   direction=-1 (lower = more active (keep min))
molecule: XJPUIXNWFNHUPM-UHFFFAOYSA-N

raw records:
   value relation
38.03375        =
38.03375        =

value kept in 21_datasets: 38.03375019908292
```

### abaumannii - step 21 preserves censoring qualifier

```
pair: CHEMBL1640427 / MIC   direction=-1
molecule: FFUICBYFZGGTEH-UHFFFAOYSA-N

raw source record   : value=696.1932371353824  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=696.1932371353824  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### abaumannii - step 22 conservative labelling (DR)

```
pair: CHEMBL5529171 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey    value relation  bin                                               reason
WFLYZIDAKLWLDY-UHFFFAOYSA-N 0.485778        =    1 '=0.4857783529393476' proves value <= 10.0 -> active
VZZYKODXNADEKN-UHFFFAOYSA-N 0.545221        =    1 '=0.5452205962532441' proves value <= 10.0 -> active
KUHQOLZBMIILJD-UHFFFAOYSA-N 1.090441        =    1 '=1.0904411925064883' proves value <= 10.0 -> active
BQJBWTDOVXRZHB-UHFFFAOYSA-N 1.161505        =    1 '=1.1615049387189995' proves value <= 10.0 -> active
QGPKADBNRMWEQR-UHFFFAOYSA-N 1.366897        =    1 '=1.3668970343801945' proves value <= 10.0 -> active
QKRZPYOEBAIZHK-UHFFFAOYSA-N 1.372480        =    1 '=1.3724801264877686' proves value <= 10.0 -> active
IFTNDGVPSSVEKH-UHFFFAOYSA-N 1.978909        =    1  '=1.978908790114955' proves value <= 10.0 -> active
BBSIBLVTSWYWSL-UHFFFAOYSA-N 2.095228        =    1  '=2.095228117961343' proves value <= 10.0 -> active
```

### abaumannii - step 22 conservative labelling (SP)

```
pair: CHEMBL5153461 / ACTIVITY   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
WEKWOFNMADUEEG-UHFFFAOYSA-N   70.2        =    1 '=70.2' proves value >= 50.0 -> active
KADOFWZNMWHQLM-UHFFFAOYSA-N   69.9        =    1 '=69.9' proves value >= 50.0 -> active
ZPJNNYJNMMCUJF-UHFFFAOYSA-N   63.1        =    1 '=63.1' proves value >= 50.0 -> active
RKAIGHYWXYTFFA-UHFFFAOYSA-N   59.3        =    1 '=59.3' proves value >= 50.0 -> active
NAAGYSIXSICQHU-UHFFFAOYSA-N   53.1        =    1 '=53.1' proves value >= 50.0 -> active
ZVZGPIABJGVAGE-UHFFFAOYSA-N   49.5        =    0   value 49.5 beyond cutoff -> inactive
OMUGTCDWYXDLAZ-UHFFFAOYSA-N   46.6        =    0   value 46.6 beyond cutoff -> inactive
HZCGMILYNXSJJT-UHFFFAOYSA-N   44.9        =    0   value 44.9 beyond cutoff -> inactive
```

### calbicans - step 21 dedup keeps most-active record

```
pair: CHEMBL1614164 / EC50   direction=-1 (lower = more active (keep min))
molecule: RNULTJAJPNAXDW-UHFFFAOYSA-M

raw records:
 value relation
 0.227        =
 0.290        =

value kept in 21_datasets: 0.227
```

### calbicans - step 21 preserves censoring qualifier

```
pair: CHEMBL1614164 / EC50   direction=-1
molecule: MNULQEHRPTZRPB-UHFFFAOYSA-N

raw source record   : value=38.0  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=38.0  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### calbicans - step 22 conservative labelling (DR)

```
pair: CHEMBL1614164 / EC50   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
FPCBQKODIUMKDZ-UHFFFAOYSA-N   0.03        <    1 '<0.03' proves value <= 10.0 -> active
RRUFFISIXQTINO-UHFFFAOYSA-M   0.03        <    1 '<0.03' proves value <= 10.0 -> active
PMUWJCVZYFHUEV-UHFFFAOYSA-M   0.03        <    1 '<0.03' proves value <= 10.0 -> active
QMTPMEYJGNJHNG-UHFFFAOYSA-N   0.03        <    1 '<0.03' proves value <= 10.0 -> active
QWLHVUDCVJRNSA-PJQLUOCWSA-N   0.03        <    1 '<0.03' proves value <= 10.0 -> active
ZYSHYUDNGIPRBJ-UHFFFAOYSA-N   0.03        <    1 '<0.03' proves value <= 10.0 -> active
QSIWCXBQJSETIC-UHFFFAOYSA-N   0.03        <    1 '<0.03' proves value <= 10.0 -> active
XIJOBJBWOPBBNG-UHFFFAOYSA-N   0.03        <    1 '<0.03' proves value <= 10.0 -> active
```

### calbicans - step 22 conservative labelling (SP)

```
pair: CHEMBL3077839 / IZ   category=SP   unit=mm   direction=1   cutoff=20.0
rule: active iff relation in {=,>} and value >= 20.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
XMAYWYJOQHXEEK-UHFFFAOYSA-N   33.0        =    1 '=33.0' proves value >= 20.0 -> active
GKIHITDQPNQOBC-UHFFFAOYSA-N   28.0        =    1 '=28.0' proves value >= 20.0 -> active
HSLMTBJBOJRELZ-UHFFFAOYSA-N   27.0        =    1 '=27.0' proves value >= 20.0 -> active
OPDYGBGGVQXNQR-UHFFFAOYSA-N   26.0        =    1 '=26.0' proves value >= 20.0 -> active
VSASCJOIQJUJKN-UHFFFAOYSA-N   26.0        =    1 '=26.0' proves value >= 20.0 -> active
QMNRBBCYLRFDIG-UHFFFAOYSA-N   25.0        =    1 '=25.0' proves value >= 20.0 -> active
GIEWFQIZOJJLQV-UHFFFAOYSA-N   24.0        =    1 '=24.0' proves value >= 20.0 -> active
DLQHTIFWTPVZGM-UHFFFAOYSA-N   24.0        =    1 '=24.0' proves value >= 20.0 -> active
```

### campylobacter - step 21 preserves censoring qualifier

```
pair: CHEMBL3117430 / MIC   direction=-1
molecule: YESFMCKAHOVAGH-VCPQQGNGSA-N

raw source record   : value=64.31836044781015  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=64.31836044781015  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### campylobacter - step 22 conservative labelling (DR)

```
pair: CHEMBL3117430 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey    value relation  bin                                               reason
PLCAYWJCESQEBF-ZAUGLNQOSA-N 0.253662        =    1 '=0.2536616052721028' proves value <= 10.0 -> active
FSWZYRGQRGYQRR-IASSRSGUSA-N 0.255938        =    1 '=0.2559382472197108' proves value <= 10.0 -> active
JKOSDVHRHMBPPA-YCVAMNHMSA-N 0.258256        =    1  '=0.258256125512477' proves value <= 10.0 -> active
SWZWUNILWPYIQU-LVDBLVHASA-N 0.482423        =    1  '=0.482422503154301' proves value <= 10.0 -> active
AVTXPUXBULDSGU-UORHQLCBSA-N 0.508202        =    1 '=0.5082019229578916' proves value <= 10.0 -> active
UIADEIUIJBZENF-XMKRHESVSA-N 0.512889        =    1 '=0.5128893027072664' proves value <= 10.0 -> active
GKOHGRAOENXRJJ-HCXNLQMVSA-N 0.958891        =    1 '=0.9588912266984824' proves value <= 10.0 -> active
GTSTTWKWLIZPCC-AOWXSCFISA-N 0.966023        =    1  '=0.966023229143187' proves value <= 10.0 -> active
```

### campylobacter - step 22 conservative labelling (SP)

```
pair: CHEMBL1266697 / ACTIVITY   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
MHWLWQUZZRMNGJ-UHFFFAOYSA-N   83.3        =    1 '=83.3' proves value >= 50.0 -> active
MYSWGUAQZAJSOK-UHFFFAOYSA-N   83.3        =    1 '=83.3' proves value >= 50.0 -> active
OGJPXUAPXNRGGI-UHFFFAOYSA-N   83.3        =    1 '=83.3' proves value >= 50.0 -> active
GSDSWSVVBLHKDQ-UHFFFAOYSA-N   83.3        =    1 '=83.3' proves value >= 50.0 -> active
SUIQUYDRLGGZOL-RCWTXCDDSA-N   72.2        =    1 '=72.2' proves value >= 50.0 -> active
OFVLGDICTFRJMM-WESIUVDSSA-N   66.7        =    1 '=66.7' proves value >= 50.0 -> active
ULGZDMOVFRHVEP-RWJQBGPGSA-N   50.0        =    1 '=50.0' proves value >= 50.0 -> active
KDLRVYVGXIQJDK-AWPVFWJPSA-N   50.0        =    1 '=50.0' proves value >= 50.0 -> active
```

### ecoli - step 21 dedup keeps most-active record

```
pair: CHEMBL1614169 / EC50   direction=-1 (lower = more active (keep min))
molecule: LWYXFDXUMVEZKS-ZVFOLQIPSA-N

raw records:
 value relation
 110.0        >
 110.0        >

value kept in 21_datasets: 110.0
```

### ecoli - step 21 preserves censoring qualifier

```
pair: CHEMBL1614169 / EC50   direction=-1
molecule: ATNOETPPUJCYBV-QDEBKDIKSA-N

raw source record   : value=110.0  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=110.0  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### ecoli - step 22 conservative labelling (DR)

```
pair: CHEMBL1614169 / EC50   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                  reason
HXOOESROJHMMSY-FNORWQNLSA-N  0.232        =    1 '=0.232' proves value <= 10.0 -> active
RTXOFQZKPXMALH-GHXIOONMSA-N  0.521        =    1 '=0.521' proves value <= 10.0 -> active
RRXFFFRYNLDYQX-BEQNREIPSA-N  4.585        =    1 '=4.585' proves value <= 10.0 -> active
JRNIHERUNQWMMW-WBYAVNBMSA-N  4.674        =    1 '=4.674' proves value <= 10.0 -> active
WTJXVDPDEQKTCV-GBTZMIKPSA-N  6.493        =    1 '=6.493' proves value <= 10.0 -> active
WNPLITIXMFABAF-WOYDFFIASA-N  7.378        =    1 '=7.378' proves value <= 10.0 -> active
KTUXFHGUHYRDNF-UHFFFAOYSA-N  8.681        =    1 '=8.681' proves value <= 10.0 -> active
QPERNWXRUSFVPU-UHFFFAOYSA-M  9.414        =    1 '=9.414' proves value <= 10.0 -> active
```

### ecoli - step 22 conservative labelling (SP)

```
pair: CHEMBL910461 / IZ   category=SP   unit=mm   direction=1   cutoff=20.0
rule: active iff relation in {=,>} and value >= 20.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
WIIZWVCIJKGZOK-RKDXNWHRSA-N   24.0        =    1 '=24.0' proves value >= 20.0 -> active
OFVLGDICTFRJMM-WESIUVDSSA-N   21.0        =    1 '=21.0' proves value >= 20.0 -> active
SBUJHOSQTJFQJX-NOAMYHISSA-N   21.0        =    1 '=21.0' proves value >= 20.0 -> active
GJXWDTUCERCKIX-UHFFFAOYSA-N   20.0        =    1 '=20.0' proves value >= 20.0 -> active
VFLFYDXGZMGUAF-UHFFFAOYSA-N   19.0        =    0   value 19.0 beyond cutoff -> inactive
PGBHMTALBVVCIT-VCIWKGPPSA-N   15.0        =    0   value 15.0 beyond cutoff -> inactive
UCSJYZPVAKXKNQ-HZYVHMACSA-N   15.0        =    0   value 15.0 beyond cutoff -> inactive
YETRAIFXPZXEDB-UHFFFAOYSA-N   14.0        =    0   value 14.0 beyond cutoff -> inactive
```

### efaecium - step 21 dedup keeps most-active record

```
pair: CHEMBL840935 / MIC   direction=-1 (lower = more active (keep min))
molecule: JNISLSXAMCNFHI-DNVCBOLYSA-N

raw records:
     value relation
233.757933        =
233.757933        >
233.757933        >
233.757933        >
233.757933        >
233.757933        >
233.757933        >
233.757933        >
233.757933        >
233.757933        >

value kept in 21_datasets: 233.75793270328265
```

### efaecium - step 21 preserves censoring qualifier

```
pair: CHEMBL813256 / MIC   direction=-1
molecule: ZWOXSWZFCWWARU-WKLXBBMDSA-N

raw source record   : value=126.69579350169356  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=126.69579350169356  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### efaecium - step 22 conservative labelling (DR)

```
pair: CHEMBL813256 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey    value relation  bin                                               reason
CKPUCCJVTDZOTQ-CBLAECCRSA-N 0.262009        =    1 '=0.2620094658779832' proves value <= 10.0 -> active
NHTGXFKVXVYFKS-QWQAMPRASA-N 0.530699        =    1 '=0.5306993662388169' proves value <= 10.0 -> active
YYDGIGWSUIKDEC-BNLPPSSXSA-N 0.530699        =    1 '=0.5306993662388169' proves value <= 10.0 -> active
KVLMFFYXICQREI-KEIFWEBJSA-N 0.539867        =    1 '=0.5398670199556446' proves value <= 10.0 -> active
OWTBXVFHEMDVHO-CCGYTOTOSA-N 0.550605        =    1 '=0.5506050047792515' proves value <= 10.0 -> active
BXCAOHIJTXVWHN-FYXRSSPKSA-N 0.688784        =    1 '=0.6887839794356655' proves value <= 10.0 -> active
FERSDKADYZRIAA-CQGKBTLCSA-N 0.859480        =    1 '=0.8594801520248494' proves value <= 10.0 -> active
APLXZEVIXWPRCD-CCGYTOTOSA-N 1.102405        =    1  '=1.102404564836822' proves value <= 10.0 -> active
```

### efaecium - step 22 conservative labelling (SP)

```
pair: CHEMBL4666314 / INHIBITION   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
YWJJMGOABBHKMZ-UHFFFAOYSA-N   97.0        =    1 '=97.0' proves value >= 50.0 -> active
TYUGIMHBLIPLKX-SNVBAGLBSA-N   79.0        =    1 '=79.0' proves value >= 50.0 -> active
XTYGWIYPRQZJMA-HNNXBMFYSA-N   73.0        =    1 '=73.0' proves value >= 50.0 -> active
CEDOVDFAPKPYNH-AWEZNQCLSA-N   67.0        =    1 '=67.0' proves value >= 50.0 -> active
CEDOVDFAPKPYNH-CQSZACIVSA-N   59.0        =    1 '=59.0' proves value >= 50.0 -> active
HYVHZJDRXRFBAL-UHFFFAOYSA-N   54.0        =    1 '=54.0' proves value >= 50.0 -> active
JJNJJFXIHCARFH-UHFFFAOYSA-N   53.0        =    1 '=53.0' proves value >= 50.0 -> active
XTYGWIYPRQZJMA-OAHLLOKOSA-N   52.0        =    1 '=52.0' proves value >= 50.0 -> active
```

### enterobacter - step 21 dedup keeps most-active record

```
pair: CHEMBL679038 / MIC   direction=-1 (lower = more active (keep min))
molecule: SARMVEMODGTEEG-WDEREUQCSA-N

raw records:
   value relation
0.267083        =
0.267083        =

value kept in 21_datasets: 0.2670826033075509
```

### enterobacter - step 21 preserves censoring qualifier

```
pair: CHEMBL4194366 / MIC   direction=-1
molecule: ACSXSYCLVXXHIZ-UHFFFAOYSA-N

raw source record   : value=1495.1087750036502  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=1495.1087750036502  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### enterobacter - step 22 conservative labelling (DR)

```
pair: CHEMBL857113 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey    value relation  bin                                               reason
YPSUKVYDKBNNFM-UHFFFAOYSA-N 0.035396        =    1 '=0.0353958912449442' proves value <= 10.0 -> active
OXXAFOWSYXTTSJ-UHFFFAOYSA-N 0.035517        =    1 '=0.0355165767702645' proves value <= 10.0 -> active
FRGTYWHQLROWDS-UHFFFAOYSA-N 0.040347        =    1 '=0.0403474181687113' proves value <= 10.0 -> active
PTMTWDHCSVIFST-UHFFFAOYSA-N 0.041007        =    1 '=0.0410069110314058' proves value <= 10.0 -> active
WXMYDHQMJVPNTD-UHFFFAOYSA-N 0.041169        =    1 '=0.0411689794484454' proves value <= 10.0 -> active
ARRLOGFHLGTCLM-UHFFFAOYSA-N 0.041619        =    1 '=0.0416188050408696' proves value <= 10.0 -> active
IUQKXTKLHSXZNO-UHFFFAOYSA-N 0.042642        =    1 '=0.0426421048142936' proves value <= 10.0 -> active
KUBSPTAIVIRIHC-UHFFFAOYSA-N 0.042817        =    1 '=0.0428173838578462' proves value <= 10.0 -> active
```

### enterobacter - step 22 conservative labelling (SP)

```
pair: CHEMBL4393348 / IZ   category=SP   unit=mm   direction=1   cutoff=20.0
rule: active iff relation in {=,>} and value >= 20.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
HHDZHSAKQQVWBJ-JPYJTQIMSA-N   36.0        =    1 '=36.0' proves value >= 20.0 -> active
ADYOKYBYYYORLP-DEOSSOPVSA-N   20.0        =    1 '=20.0' proves value >= 20.0 -> active
CBIQJLKZPYDEDQ-QHCPKHFHSA-N   14.0        =    0   value 14.0 beyond cutoff -> inactive
KDJBGLNOJOAKCM-GOTSBHOMSA-N   14.0        =    0   value 14.0 beyond cutoff -> inactive
JWOWBQWOUVZHIO-DEOSSOPVSA-N   14.0        =    0   value 14.0 beyond cutoff -> inactive
QRQKVQCDKFEWSW-QHCPKHFHSA-N   13.0        =    0   value 13.0 beyond cutoff -> inactive
JIEZBUYCEJCVLP-DEOSSOPVSA-N   10.0        =    0   value 10.0 beyond cutoff -> inactive
VVGZNYKBTPSOQC-QHCPKHFHSA-N   10.0        =    0   value 10.0 beyond cutoff -> inactive
```

### hpylori - step 21 dedup keeps most-active record

```
pair: CHEMBL1005166 / MIC   direction=-1 (lower = more active (keep min))
molecule: MDVWGBKYBVEXHU-UHFFFAOYSA-N

raw records:
   value relation
0.047786        =
0.047786        =

value kept in 21_datasets: 0.0477862996678852
```

### hpylori - step 21 preserves censoring qualifier

```
pair: CHEMBL1291866 / MIC   direction=-1
molecule: VTIPSNMJULYGDJ-UHFFFAOYSA-N

raw source record   : value=192.1591385121502  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=192.1591385121502  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### hpylori - step 22 conservative labelling (DR)

```
pair: CHEMBL1291866 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey     value relation  bin                                               reason
TUTMEHZVWWQNMW-UHFFFAOYSA-N  1.948120        =    1  '=1.948120064138107' proves value <= 10.0 -> active
VAOCPAMSLUNLGC-UHFFFAOYSA-N  2.921312        =    1  '=2.921311552034401' proves value <= 10.0 -> active
NZGZWRPXPMIHRW-UHFFFAOYSA-N  3.607904        =    1  '=3.607904085257549' proves value <= 10.0 -> active
DHJPGNFCNCQUMP-UHFFFAOYSA-N  5.589615        =    1  '=5.589614782381248' proves value <= 10.0 -> active
CFMUZTVJQKVMAZ-UHFFFAOYSA-N  5.972686        =    1  '=5.972686445447549' proves value <= 10.0 -> active
VGBSVSAPEVDVLS-UHFFFAOYSA-N  6.041150        =    1 '=6.0411495223232095' proves value <= 10.0 -> active
HEPOSDYQFQRNIC-UHFFFAOYSA-N  7.283729        =    1  '=7.283729269385926' proves value <= 10.0 -> active
VVLFUNBKTUJMPN-UHFFFAOYSA-N 13.496970        =    0   value 13.496969670751806 beyond cutoff -> inactive
```

### hpylori - step 22 conservative labelling (SP)

```
pair: CHEMBL2166540 / IZ   category=SP   unit=mm   direction=1   cutoff=20.0
rule: active iff relation in {=,>} and value >= 20.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
LUDIEESNHSPSDG-UHFFFAOYSA-N   24.0        =    1 '=24.0' proves value >= 20.0 -> active
KTKBBFFRASJHLH-UHFFFAOYSA-N   24.0        =    1 '=24.0' proves value >= 20.0 -> active
IKDQPJVFRIGXLA-UHFFFAOYSA-N   22.0        =    1 '=22.0' proves value >= 20.0 -> active
FQZBULWMJLBFPY-UHFFFAOYSA-N   22.0        =    1 '=22.0' proves value >= 20.0 -> active
YWNRRGHBQRUPQE-UHFFFAOYSA-N   21.0        =    1 '=21.0' proves value >= 20.0 -> active
BOCYDKUJDBGHER-UHFFFAOYSA-N   21.0        =    1 '=21.0' proves value >= 20.0 -> active
QLBDQAVJICCBNR-UHFFFAOYSA-N   21.0        =    1 '=21.0' proves value >= 20.0 -> active
OAGIORMWSLAGIO-UHFFFAOYSA-N   19.0        =    0   value 19.0 beyond cutoff -> inactive
```

### kpneumoniae - step 21 dedup keeps most-active record

```
pair: CHEMBL700126 / MIC   direction=-1 (lower = more active (keep min))
molecule: SARMVEMODGTEEG-WDEREUQCSA-N

raw records:
   value relation
0.534165        =
0.534165        =

value kept in 21_datasets: 0.5341652066151019
```

### kpneumoniae - step 21 preserves censoring qualifier

```
pair: CHEMBL5129403 / MIC   direction=-1
molecule: TYZROVQLWOKYKF-ZDUSSCGKSA-N

raw source record   : value=94.85669228785449  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=94.85669228785449  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### kpneumoniae - step 22 conservative labelling (DR)

```
pair: CHEMBL5129403 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey    value relation  bin                                               reason
DMJNNHOOLUXYBV-PQTSNVLCSA-N 0.162985        <    1 '<0.1629853704331499' proves value <= 10.0 -> active
HXAZKAMWTTVQTI-XAVRLVDKSA-N 1.091443        =    1 '=1.0914433027947497' proves value <= 10.0 -> active
IHJKKUFHLWJRNI-HFZQUGABSA-N 1.108413        =    1  '=1.108412741869515' proves value <= 10.0 -> active
XQQSXVFRSJIKNI-JOMAPVKXSA-N 1.141353        =    1 '=1.1413531654859368' proves value <= 10.0 -> active
UPEHOGVLCWTGAR-NLCGWVTRSA-N 2.251836        =    1  '=2.251836372561824' proves value <= 10.0 -> active
AZRUCCFTSHLJLH-LDMPPBJASA-N 2.358215        =    1 '=2.3582152555303093' proves value <= 10.0 -> active
VPVQSKONIULWRI-GDLQSTFJSA-N 2.435968        =    1  '=2.435967548040325' proves value <= 10.0 -> active
SLCPVXQXKOPIFE-YABOCYAYSA-N 3.823046        =    1 '=3.8230464710413785' proves value <= 10.0 -> active
```

### kpneumoniae - step 22 conservative labelling (SP)

```
pair: CHEMBL1285737 / IZ   category=SP   unit=mm   direction=1   cutoff=20.0
rule: active iff relation in {=,>} and value >= 20.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
MYSWGUAQZAJSOK-UHFFFAOYSA-N   22.0        =    1 '=22.0' proves value >= 20.0 -> active
KEAJPYMVAICCDZ-UHFFFAOYSA-N   16.0        =    0   value 16.0 beyond cutoff -> inactive
SWUHOHQXIOJTKN-UHFFFAOYSA-N   16.0        =    0   value 16.0 beyond cutoff -> inactive
BNYWSAYALSCTGS-UHFFFAOYSA-N   16.0        =    0   value 16.0 beyond cutoff -> inactive
ZOBYCNABTAKCBM-UHFFFAOYSA-N   16.0        =    0   value 16.0 beyond cutoff -> inactive
OSIZKECJWFLJER-UHFFFAOYSA-N   16.0        =    0   value 16.0 beyond cutoff -> inactive
USJRUYQKZLXDPT-UHFFFAOYSA-N   16.0        =    0   value 16.0 beyond cutoff -> inactive
CWGGTXKORNGWCY-UHFFFAOYSA-N   16.0        =    0   value 16.0 beyond cutoff -> inactive
```

### mtuberculosis - step 21 dedup keeps most-active record

```
pair: CHEMBL4649948 / PERCENTEFFECT   direction=1 (higher = more active (keep max))
molecule: JMOQCWDQWVVQKQ-UHFFFAOYSA-N

raw records:
  value relation
-10.680        =
 -2.002        =

value kept in 21_datasets: -2.002
```

### mtuberculosis - step 21 preserves censoring qualifier

```
pair: CHEMBL4649949 / IC50   direction=-1
molecule: JDYLTJYGARCQNF-UHFFFAOYSA-N

raw source record   : value=50.11872  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=50.11872  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### mtuberculosis - step 22 conservative labelling (DR)

```
pair: CHEMBL4649949 / IC50   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey   value relation  bin                                               reason
ACTIMEXECSWMLS-UHFFFAOYSA-N 0.19778        =    1            '=0.19778' proves value <= 10.0 -> active
AJQOZMQCHSYLDT-UHFFFAOYSA-N 0.21091        =    1            '=0.21091' proves value <= 10.0 -> active
YQQXDQMDDXMYAP-UHFFFAOYSA-N 0.60922        =    1            '=0.60922' proves value <= 10.0 -> active
VDPMCMAKCOQWEF-UHFFFAOYSA-N 0.70415        =    1 '=0.7041499999999999' proves value <= 10.0 -> active
GNIMELCEHUTLKM-UHFFFAOYSA-N 0.88375        =    1            '=0.88375' proves value <= 10.0 -> active
RDYBQOMUCATMBF-UHFFFAOYSA-N 0.89282        =    1 '=0.8928200000000001' proves value <= 10.0 -> active
KUOBCWMMYIOLIA-UHFFFAOYSA-N 0.92810        =    1             '=0.9281' proves value <= 10.0 -> active
LMJXZAXRPCVNLB-UHFFFAOYSA-N 0.95705        =    1            '=0.95705' proves value <= 10.0 -> active
```

### mtuberculosis - step 22 conservative labelling (SP)

```
pair: CHEMBL4649948 / PERCENTEFFECT   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                   reason
CDRBRQSPQWOTSS-UHFFFAOYSA-N 120.27        =    1 '=120.27' proves value >= 50.0 -> active
AFSHJZFNKRCNER-UHFFFAOYSA-N 118.49        =    1 '=118.49' proves value >= 50.0 -> active
LGYWZWVXJPZZOI-UHFFFAOYSA-N 117.34        =    1 '=117.34' proves value >= 50.0 -> active
CAOAXNRVWBYUSK-LKLBIAKRSA-N 117.30        =    1  '=117.3' proves value >= 50.0 -> active
XYMOGNZVKMQGSM-UHFFFAOYSA-N 116.09        =    1 '=116.09' proves value >= 50.0 -> active
SJIXMBBQRZPUMK-UHFFFAOYSA-N 115.96        =    1 '=115.96' proves value >= 50.0 -> active
FGEZLEBKHSAYSJ-UHFFFAOYSA-N 115.88        =    1 '=115.88' proves value >= 50.0 -> active
QYOJZFBQEAZNEW-UHFFFAOYSA-N 115.84        =    1 '=115.84' proves value >= 50.0 -> active
```

### ngonorrhoeae - step 21 preserves censoring qualifier

```
pair: CHEMBL753615 / MIC   direction=-1
molecule: MPYBCPFIBJZKPK-UHFFFAOYSA-N

raw source record   : value=31.60905915635421  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=31.60905915635421  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### ngonorrhoeae - step 22 conservative labelling (DR)

```
pair: CHEMBL753615 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey    value relation  bin                                               reason
JSRFRMHDGUSXFV-UHFFFAOYSA-N 0.595922        =    1 '=0.5959219866608823' proves value <= 10.0 -> active
IIHNEDBVMJJISR-UHFFFAOYSA-N 0.750203        =    1 '=0.7502033050956809' proves value <= 10.0 -> active
QQDFSPOTSJZPSS-UHFFFAOYSA-N 1.090957        =    1 '=1.0909574897504544' proves value <= 10.0 -> active
CCUMMLKYWTZUTM-UHFFFAOYSA-N 1.233073        =    1  '=1.233072990522601' proves value <= 10.0 -> active
XJGQPBMAGXHFHK-UHFFFAOYSA-N 1.276666        =    1 '=1.2766663687778474' proves value <= 10.0 -> active
NKZRDMDWAYJJNS-UHFFFAOYSA-N 1.298455        =    1 '=1.2984550981242515' proves value <= 10.0 -> active
FYBBPDOWFMLMIZ-UHFFFAOYSA-N 1.298455        =    1 '=1.2984550981242515' proves value <= 10.0 -> active
SUGUMJSHCAGFER-UHFFFAOYSA-N 1.409845        =    1  '=1.409844663314996' proves value <= 10.0 -> active
```

### ngonorrhoeae - step 22 conservative labelling (SP)

```
pair: CHEMBL5335036 / ACTIVITY   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
IILINGWUQRJNCB-LOSJGSFVSA-N   73.2        =    1 '=73.2' proves value >= 50.0 -> active
NVKAOHZJXXDLFP-SVBPBHIXSA-N   71.5        =    1 '=71.5' proves value >= 50.0 -> active
LMJKJLOIXRKFPW-FTJBHMTQSA-N   70.6        =    1 '=70.6' proves value >= 50.0 -> active
IILINGWUQRJNCB-RPBOFIJWSA-N   68.4        =    1 '=68.4' proves value >= 50.0 -> active
VNXURXHRNFTYAY-AMVUTOCUSA-N   65.1        =    1 '=65.1' proves value >= 50.0 -> active
HBCVEDVLIMIBJN-KEKNWZKVSA-N   65.1        =    1 '=65.1' proves value >= 50.0 -> active
IILINGWUQRJNCB-DQEYMECFSA-N   64.6        =    1 '=64.6' proves value >= 50.0 -> active
HLORBPIXFLEOFV-CPRJBALCSA-N   63.5        =    1 '=63.5' proves value >= 50.0 -> active
```

### paeruginosa - step 21 dedup keeps most-active record

```
pair: CHEMBL1640426 / MIC   direction=-1 (lower = more active (keep min))
molecule: XJPUIXNWFNHUPM-UHFFFAOYSA-N

raw records:
   value relation
38.03375        =
38.03375        =

value kept in 21_datasets: 38.03375019908292
```

### paeruginosa - step 21 preserves censoring qualifier

```
pair: CHEMBL1640426 / MIC   direction=-1
molecule: FFUICBYFZGGTEH-UHFFFAOYSA-N

raw source record   : value=696.1932371353824  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=696.1932371353824  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### paeruginosa - step 22 conservative labelling (DR)

```
pair: CHEMBL1640426 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey     value relation  bin                                              reason
PPWMMWQBJRDQFH-SOFGYWHQSA-N  5.859552        =    1 '=5.859552388793021' proves value <= 10.0 -> active
OVLKACDWJNGFOC-YVMONPNESA-N  6.703648        =    1 '=6.703648460674722' proves value <= 10.0 -> active
SBNIQRGKIQSRPT-SOFGYWHQSA-N  6.727890        =    1  '=6.72789046994315' proves value <= 10.0 -> active
CLWLYALISAAAPE-ALCCZGGFSA-N  7.034377        =    1  '=7.03437700040096' proves value <= 10.0 -> active
LXHTXRZCFIWDPF-UHFFFAOYSA-N 21.914116        =    0    value 21.9141158407556 beyond cutoff -> inactive
PTIXJNCZRFOXHN-UHFFFAOYSA-N 23.666013        =    0   value 23.66601338608882 beyond cutoff -> inactive
YRBPUSAZSOXBMT-UHFFFAOYSA-N 23.666013        =    0   value 23.66601338608882 beyond cutoff -> inactive
FYEYHDOOAYIROF-UHFFFAOYSA-N 23.736020        =    0  value 23.736019855180608 beyond cutoff -> inactive
```

### paeruginosa - step 22 conservative labelling (SP)

```
pair: CHEMBL3411950 / EFFICACY   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                  reason
JUJLIDFXDQHGLM-UHFFFAOYSA-N  103.0        =    1 '=103.0' proves value >= 50.0 -> active
SCGRZAGHUQOFNO-UHFFFAOYSA-N  101.0        =    1 '=101.0' proves value >= 50.0 -> active
JHJXLNOTXCUSJH-UHFFFAOYSA-N  100.0        =    1 '=100.0' proves value >= 50.0 -> active
RBYAEXHPLASZEL-UHFFFAOYSA-N   99.0        =    1  '=99.0' proves value >= 50.0 -> active
XKYQAWFNOGSGSZ-UHFFFAOYSA-N   97.0        =    1  '=97.0' proves value >= 50.0 -> active
NJIKTUINATYSPX-UHFFFAOYSA-N   96.0        =    1  '=96.0' proves value >= 50.0 -> active
QFAZPNYVHPVYMY-UHFFFAOYSA-N   95.0        =    1  '=95.0' proves value >= 50.0 -> active
XOTDJPZBRPPQKF-UHFFFAOYSA-N   95.0        =    1  '=95.0' proves value >= 50.0 -> active
```

### pfalciparum - step 21 dedup keeps most-active record

```
pair: CHEMBL4888485 / INHIBITION   direction=1 (higher = more active (keep max))
molecule: BITOFCODWBHKHG-UHFFFAOYSA-N

raw records:
 value relation
   1.0        =
  -9.0        =
  -4.0        =

value kept in 21_datasets: 1.0
```

### pfalciparum - step 21 preserves censoring qualifier

```
pair: CHEMBL1040691 / EC50   direction=-1
molecule: SXPLQBQKDMNRHS-UHFFFAOYSA-N

raw source record   : value=1.25  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=1.25  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### pfalciparum - step 22 conservative labelling (DR)

```
pair: CHEMBL1040691 / EC50   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 1082

                   inchikey    value relation  bin                                               reason
NSFWWJIQIKBZMJ-YBIYSQHJSA-N 0.000572        <    1 '<0.0005719999999999' proves value <= 10.0 -> active
FYJWEIIPLUAXTN-UHFFFAOYSA-N 0.000592        =    1           '=0.000592' proves value <= 10.0 -> active
MOUAPRKJJUXEIE-UHFFFAOYSA-N 0.000986        =    1           '=0.000986' proves value <= 10.0 -> active
PPMITYXUFZYPTA-UHFFFAOYSA-N 0.001231        =    1           '=0.001231' proves value <= 10.0 -> active
WHTVZRBIWZFKQO-UHFFFAOYSA-N 0.002426        =    1           '=0.002426' proves value <= 10.0 -> active
BAQQTXLPROBUPE-UHFFFAOYSA-N 0.002428        =    1           '=0.002428' proves value <= 10.0 -> active
MBAIJEGYSUBFHC-UHFFFAOYSA-N 0.002431        =    1           '=0.002431' proves value <= 10.0 -> active
BOMPKFLFSCDGLF-HFXFMHOMSA-N 0.002435        =    1           '=0.002435' proves value <= 10.0 -> active
```

### pfalciparum - step 22 conservative labelling (SP)

```
pair: CHEMBL4888485 / INHIBITION   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                  reason
NEVNBUVUQAQGTC-UHFFFAOYSA-N  105.0        =    1 '=105.0' proves value >= 50.0 -> active
AGUDRLRAJMJXQN-UHFFFAOYSA-N  101.0        =    1 '=101.0' proves value >= 50.0 -> active
OJZZOAHEPFEYPB-UHFFFAOYSA-N  101.0        =    1 '=101.0' proves value >= 50.0 -> active
CQTGJKYBUCHLFO-UHFFFAOYSA-N  100.0        =    1 '=100.0' proves value >= 50.0 -> active
SNJZJAFGTKFQKP-UHFFFAOYSA-N  100.0        =    1 '=100.0' proves value >= 50.0 -> active
CFJWEEFKGCETGA-MGPLVRAMSA-N  100.0        =    1 '=100.0' proves value >= 50.0 -> active
LEWOIHLFDUTZAO-UHFFFAOYSA-N   99.0        =    1  '=99.0' proves value >= 50.0 -> active
BNEUIECAJZZTBO-UHFFFAOYSA-N   99.0        =    1  '=99.0' proves value >= 50.0 -> active
```

### saureus - step 21 dedup keeps most-active record

```
pair: CHEMBL1640423 / MIC   direction=-1 (lower = more active (keep min))
molecule: XJPUIXNWFNHUPM-UHFFFAOYSA-N

raw records:
   value relation
9.508438        >
0.071313        <

value kept in 21_datasets: 0.0713132816232804
```

### saureus - step 21 preserves censoring qualifier

```
pair: CHEMBL1640423 / MIC   direction=-1
molecule: JPONBNWNHDAQOF-UHFFFAOYSA-N

raw source record   : value=176.56341385360685  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=176.56341385360685  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### saureus - step 22 conservative labelling (DR)

```
pair: CHEMBL1640423 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 1

                   inchikey    value relation  bin                                               reason
PTIXJNCZRFOXHN-UHFFFAOYSA-N 0.044374        <    1 '<0.0443737750989165' proves value <= 10.0 -> active
YRBPUSAZSOXBMT-UHFFFAOYSA-N 0.044374        <    1 '<0.0443737750989165' proves value <= 10.0 -> active
FYEYHDOOAYIROF-UHFFFAOYSA-N 0.044505        <    1 '<0.0445050372284636' proves value <= 10.0 -> active
KWWKOJDMVGWSGT-UHFFFAOYSA-N 0.046319        <    1 '<0.0463193129919497' proves value <= 10.0 -> active
ZTSOCMPHOSIHFE-UHFFFAOYSA-N 0.047020        <    1 '<0.0470200273970026' proves value <= 10.0 -> active
JQTVADHPYHXVJZ-UHFFFAOYSA-N 0.047639        <    1 '<0.0476393132316604' proves value <= 10.0 -> active
OVHNRTXDZSQCJF-UHFFFAOYSA-N 0.047689        <    1 '<0.0476890672813257' proves value <= 10.0 -> active
FWMPROYARROAKG-UHFFFAOYSA-N 0.047808        <    1 '<0.0478076983143005' proves value <= 10.0 -> active
```

### saureus - step 22 conservative labelling (SP)

```
pair: CHEMBL4379067 / ACTIVITY   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 5

                   inchikey  value relation  bin                                  reason
KSFAZTRLNGVKCI-GCPFILQVSA-N  110.4        =    1 '=110.4' proves value >= 50.0 -> active
HAKYAXOGASYZIE-MUOSHHONSA-N  108.6        =    1 '=108.6' proves value >= 50.0 -> active
NFSCQCDEWHYJGR-OBNUCGLBSA-N  108.3        =    1 '=108.3' proves value >= 50.0 -> active
XCVUAIHWQQCXSE-QZPDNNCWSA-N  107.6        =    1 '=107.6' proves value >= 50.0 -> active
ZSZJKDJQAIFACO-GCPFILQVSA-N  107.0        =    1 '=107.0' proves value >= 50.0 -> active
RDAZGWDWLZHKSO-XXQCFZTBSA-N  106.7        =    1 '=106.7' proves value >= 50.0 -> active
XULRGMMZGKTEMZ-DZOYZYGOSA-N  105.9        =    1 '=105.9' proves value >= 50.0 -> active
PZTILDNQKFDUDZ-DCLMVBETSA-N  104.8        =    1 '=104.8' proves value >= 50.0 -> active
```

### smansoni - step 21 preserves censoring qualifier

```
pair: CHEMBL4651365 / EC50   direction=-1
molecule: FNJYTAKMPPRBPA-UHFFFAOYSA-N

raw source record   : value=25.0  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=25.0  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### smansoni - step 22 conservative labelling (DR)

```
pair: CHEMBL4651365 / EC50   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey   value relation  bin                                               reason
WDXWJQHSWIRCAG-UHFFFAOYSA-N 0.00267        <    1            '<0.00267' proves value <= 10.0 -> active
HMJOUKNDLMKCFL-UHFFFAOYSA-N 0.00410        =    1 '=0.0040999999999999' proves value <= 10.0 -> active
FFQPWVRLUPXYEI-UHFFFAOYSA-N 0.00430        =    1             '=0.0043' proves value <= 10.0 -> active
XGYPQHLGBHDJBR-UHFFFAOYSA-N 0.00496        =    1            '=0.00496' proves value <= 10.0 -> active
JYBCFMXGGAMLEF-UHFFFAOYSA-N 0.00518        =    1            '=0.00518' proves value <= 10.0 -> active
XPRZITWAQRBCAV-UHFFFAOYSA-N 0.00519        =    1            '=0.00519' proves value <= 10.0 -> active
ZRBZDGAFGGBUDP-UHFFFAOYSA-N 0.00540        =    1             '=0.0054' proves value <= 10.0 -> active
AAVPNSAQPOLIIG-UHFFFAOYSA-N 0.00821        =    1            '=0.00821' proves value <= 10.0 -> active
```

### smansoni - step 22 conservative labelling (SP)

```
pair: CHEMBL4649967 / PERCENTEFFECT   category=SP   unit=%   direction=1   cutoff=50.0
rule: active iff relation in {=,>} and value >= 50.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                   reason
NSHVJJLMZBEBQL-UHFFFAOYSA-N 106.57        =    1 '=106.57' proves value >= 50.0 -> active
DKHCPVUGZSMKND-UHFFFAOYSA-N 102.73        =    1 '=102.73' proves value >= 50.0 -> active
ONJKCVKAFMMXNI-UHFFFAOYSA-N  98.89        =    1  '=98.89' proves value >= 50.0 -> active
YTDYVORNAVTVDH-UHFFFAOYSA-N  94.46        =    1  '=94.46' proves value >= 50.0 -> active
JOLXRQOMZBFWKX-UHFFFAOYSA-N  89.81        =    1  '=89.81' proves value >= 50.0 -> active
IVKSWVAWFVKWGX-UHFFFAOYSA-N  79.00        =    1   '=79.0' proves value >= 50.0 -> active
HGVDMFLUYYBPIR-UHFFFAOYSA-N  75.77        =    1  '=75.77' proves value >= 50.0 -> active
ZYRIWNXZVXIJOZ-UHFFFAOYSA-N  75.74        =    1  '=75.74' proves value >= 50.0 -> active
```

### spneumoniae - step 21 dedup keeps most-active record

```
pair: CHEMBL808852 / MIC   direction=-1 (lower = more active (keep min))
molecule: SARMVEMODGTEEG-WDEREUQCSA-N

raw records:
   value relation
0.034721        =
0.267083        =

value kept in 21_datasets: 0.0347207384299816
```

### spneumoniae - step 21 preserves censoring qualifier

```
pair: CHEMBL812245 / MIC   direction=-1
molecule: ZKQACYFUDYDDRW-GMEXIMJCSA-N

raw source record   : value=175.84035893413267  relation='>'  (wrong-side censored non-inhibitor)
preserved in output : value=175.84035893413267  relation='>'  (qualifier kept; conservative binarisation labels it inactive unless it proves active)
```

### spneumoniae - step 22 conservative labelling (DR)

```
pair: CHEMBL812245 / MIC   category=DR   unit=umol.L-1   direction=-1   cutoff=10.0
rule: active iff relation in {=,<} and value <= 10.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey    value relation  bin                                               reason
CKPUCCJVTDZOTQ-CBLAECCRSA-N 0.062882        =    1 '=0.0628822718107159' proves value <= 10.0 -> active
OWTBXVFHEMDVHO-CCGYTOTOSA-N 0.132145        =    1 '=0.1321452011470203' proves value <= 10.0 -> active
BXCAOHIJTXVWHN-FYXRSSPKSA-N 0.165308        =    1 '=0.1653081550645597' proves value <= 10.0 -> active
FERSDKADYZRIAA-CQGKBTLCSA-N 0.206275        =    1 '=0.2062752364859638' proves value <= 10.0 -> active
NHTGXFKVXVYFKS-QWQAMPRASA-N 0.265350        =    1 '=0.2653496831194084' proves value <= 10.0 -> active
YYDGIGWSUIKDEC-BNLPPSSXSA-N 0.265350        =    1 '=0.2653496831194084' proves value <= 10.0 -> active
KVLMFFYXICQREI-KEIFWEBJSA-N 0.269934        =    1 '=0.2699335099778223' proves value <= 10.0 -> active
WBPYTXDJUQJLPQ-VMXQISHHSA-N 0.272892        =    1 '=0.2728923974361213' proves value <= 10.0 -> active
```

### spneumoniae - step 22 conservative labelling (SP)

```
pair: CHEMBL4131352 / IZ   category=SP   unit=mm   direction=1   cutoff=20.0
rule: active iff relation in {=,>} and value >= 20.0
inactives whose value is on the active side but qualifier does not prove a hit (kept inactive, conservative): 0

                   inchikey  value relation  bin                                 reason
HMQYLLJHEKTZQY-ZMOGYAJESA-N   24.3        =    1 '=24.3' proves value >= 20.0 -> active
AVKUERGKIZMTKX-NJBDSQKTSA-N   23.8        =    1 '=23.8' proves value >= 20.0 -> active
REYOZFFLPDFTTB-DHRITJCHSA-N   23.1        =    1 '=23.1' proves value >= 20.0 -> active
TVKBSKUPCJFLTK-ZMOGYAJESA-N   22.5        =    1 '=22.5' proves value >= 20.0 -> active
CPXGRGAQHDSFKC-FSJBWODESA-N   22.3        =    1 '=22.3' proves value >= 20.0 -> active
NIAFKZVLNZNVLO-QPJQQBGISA-N   22.1        =    1 '=22.1' proves value >= 20.0 -> active
LCZXJRFYXHLVAO-BHGWPJFGSA-N   21.3        =    1 '=21.3' proves value >= 20.0 -> active
ISBRVEVNKOXBII-FSJBWODESA-N   21.3        =    1 '=21.3' proves value >= 20.0 -> active
```
