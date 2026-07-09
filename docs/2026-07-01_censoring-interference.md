# Censoring interference near activity cutoffs (steps 22/23)

Generated: 2026-07-01  
Tool: `tools/evaluate_censoring_interference.py`  
Cutoffs: low/middle/high from `config/expert_cutoffs.csv` (ORGANISM triple per activity_type/unit, as step 19; uniform across pathogens).

Metrics are over molecules whose FINAL binary label rests on a censored
('>'/'<') record that survived step-22 dedup. **Ambiguous** = the cutoff
falls inside the qualifier's true-value range (label genuinely unknown);
this is an inherent property of the data, independent of labelling policy.
**Face change** = a bare-threshold rewrite would give a different label
than the raw stated value (the interference the OLD binarisation caused).
The pipeline now labels conservatively (prove-the-hit), so ambiguous
records are all inactive - see the verification section below.

## Summary (per pathogen x category x level)

| Pathogen | Cat | Level | Datasets | Molecules | Censored | Value-changed | Ambiguous | Ambig % | Face-change | Face % |
|---|---|---|---|---|---|---|---|---|---|---|
| abaumannii | DR | low | 813 | 12372 | 4354 | 4180 | 93 | 0.752 | 64 | 0.517 |
| abaumannii | SP | low | 41 | 521 | 59 | 45 | 1 | 0.192 | 0 | 0.0 |
| calbicans | DR | low | 2231 | 40965 | 11129 | 10422 | 436 | 1.064 | 176 | 0.43 |
| calbicans | SP | low | 588 | 8422 | 509 | 423 | 5 | 0.059 | 20 | 0.237 |
| campylobacter | DR | low | 42 | 566 | 79 | 66 | 11 | 1.943 | 5 | 0.883 |
| campylobacter | SP | low | 4 | 56 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| ecoli | DR | low | 4711 | 78111 | 23855 | 21990 | 755 | 0.967 | 287 | 0.367 |
| ecoli | SP | low | 860 | 12996 | 792 | 650 | 75 | 0.577 | 0 | 0.0 |
| efaecium | DR | low | 629 | 10281 | 2501 | 2352 | 175 | 1.702 | 121 | 1.177 |
| efaecium | SP | low | 21 | 203 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| enterobacter | DR | low | 422 | 6228 | 1828 | 1680 | 84 | 1.349 | 30 | 0.482 |
| enterobacter | SP | low | 51 | 465 | 23 | 19 | 0 | 0.0 | 0 | 0.0 |
| hpylori | DR | low | 92 | 1361 | 289 | 285 | 21 | 1.543 | 1 | 0.073 |
| hpylori | SP | low | 37 | 456 | 28 | 1 | 0 | 0.0 | 0 | 0.0 |
| kpneumoniae | DR | low | 1563 | 24831 | 7400 | 7047 | 142 | 0.572 | 72 | 0.29 |
| kpneumoniae | SP | low | 291 | 4405 | 234 | 192 | 6 | 0.136 | 0 | 0.0 |
| mtuberculosis | DR | low | 2840 | 64415 | 15450 | 13789 | 1457 | 2.262 | 706 | 1.096 |
| mtuberculosis | SP | low | 410 | 234387 | 418 | 248 | 84 | 0.036 | 68 | 0.029 |
| ngonorrhoeae | DR | low | 46 | 825 | 212 | 149 | 10 | 1.212 | 8 | 0.97 |
| ngonorrhoeae | SP | low | 4 | 50 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| paeruginosa | DR | low | 2921 | 46265 | 13372 | 13090 | 424 | 0.916 | 109 | 0.236 |
| paeruginosa | SP | low | 597 | 9404 | 766 | 614 | 43 | 0.457 | 0 | 0.0 |
| pfalciparum | DR | low | 2710 | 75812 | 10741 | 10410 | 7374 | 9.727 | 538 | 0.71 |
| pfalciparum | SP | low | 155 | 405323 | 239 | 71 | 25 | 0.006 | 5 | 0.001 |
| saureus | DR | low | 7709 | 135588 | 27882 | 25274 | 1557 | 1.148 | 589 | 0.434 |
| saureus | SP | low | 1043 | 15722 | 1098 | 711 | 309 | 1.965 | 0 | 0.0 |
| smansoni | DR | low | 37 | 980 | 221 | 198 | 10 | 1.02 | 9 | 0.918 |
| smansoni | SP | low | 47 | 1936 | 15 | 8 | 3 | 0.155 | 0 | 0.0 |
| spneumoniae | DR | low | 1375 | 22589 | 4502 | 2693 | 269 | 1.191 | 81 | 0.359 |
| spneumoniae | SP | low | 92 | 832 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| abaumannii | DR | middle | 813 | 12372 | 4354 | 4180 | 49 | 0.396 | 35 | 0.283 |
| abaumannii | SP | middle | 41 | 521 | 59 | 45 | 1 | 0.192 | 1 | 0.192 |
| calbicans | DR | middle | 2231 | 40965 | 11129 | 10422 | 267 | 0.652 | 133 | 0.325 |
| calbicans | SP | middle | 588 | 8422 | 509 | 423 | 69 | 0.819 | 0 | 0.0 |
| campylobacter | DR | middle | 42 | 566 | 79 | 66 | 11 | 1.943 | 4 | 0.707 |
| campylobacter | SP | middle | 4 | 56 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| ecoli | DR | middle | 4711 | 78111 | 23855 | 21990 | 442 | 0.566 | 121 | 0.155 |
| ecoli | SP | middle | 860 | 12996 | 792 | 650 | 115 | 0.885 | 0 | 0.0 |
| efaecium | DR | middle | 629 | 10281 | 2501 | 2352 | 50 | 0.486 | 40 | 0.389 |
| efaecium | SP | middle | 21 | 203 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| enterobacter | DR | middle | 422 | 6228 | 1828 | 1680 | 71 | 1.14 | 16 | 0.257 |
| enterobacter | SP | middle | 51 | 465 | 23 | 19 | 3 | 0.645 | 0 | 0.0 |
| hpylori | DR | middle | 92 | 1361 | 289 | 285 | 21 | 1.543 | 1 | 0.073 |
| hpylori | SP | middle | 37 | 456 | 28 | 1 | 0 | 0.0 | 0 | 0.0 |
| kpneumoniae | DR | middle | 1563 | 24831 | 7400 | 7047 | 80 | 0.322 | 31 | 0.125 |
| kpneumoniae | SP | middle | 291 | 4405 | 234 | 192 | 28 | 0.636 | 0 | 0.0 |
| mtuberculosis | DR | middle | 2840 | 64415 | 15450 | 13789 | 790 | 1.226 | 314 | 0.487 |
| mtuberculosis | SP | middle | 410 | 234387 | 418 | 248 | 84 | 0.036 | 68 | 0.029 |
| ngonorrhoeae | DR | middle | 46 | 825 | 212 | 149 | 6 | 0.727 | 4 | 0.485 |
| ngonorrhoeae | SP | middle | 4 | 50 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| paeruginosa | DR | middle | 2921 | 46265 | 13372 | 13090 | 221 | 0.478 | 40 | 0.086 |
| paeruginosa | SP | middle | 597 | 9404 | 766 | 614 | 95 | 1.01 | 0 | 0.0 |
| pfalciparum | DR | middle | 2710 | 75812 | 10741 | 10410 | 3592 | 4.738 | 2983 | 3.935 |
| pfalciparum | SP | middle | 155 | 405323 | 239 | 71 | 162 | 0.04 | 0 | 0.0 |
| saureus | DR | middle | 7709 | 135588 | 27882 | 25274 | 900 | 0.664 | 141 | 0.104 |
| saureus | SP | middle | 1043 | 15722 | 1098 | 711 | 429 | 2.729 | 0 | 0.0 |
| smansoni | DR | middle | 37 | 980 | 221 | 198 | 0 | 0.0 | 9 | 0.918 |
| smansoni | SP | middle | 47 | 1936 | 15 | 8 | 8 | 0.413 | 0 | 0.0 |
| spneumoniae | DR | middle | 1375 | 22589 | 4502 | 2693 | 158 | 0.699 | 76 | 0.336 |
| spneumoniae | SP | middle | 92 | 832 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| abaumannii | DR | high | 813 | 12372 | 4354 | 4180 | 34 | 0.275 | 11 | 0.089 |
| abaumannii | SP | high | 41 | 521 | 59 | 45 | 15 | 2.879 | 1 | 0.192 |
| calbicans | DR | high | 2231 | 40965 | 11129 | 10422 | 250 | 0.61 | 116 | 0.283 |
| calbicans | SP | high | 588 | 8422 | 509 | 423 | 86 | 1.021 | 0 | 0.0 |
| campylobacter | DR | high | 42 | 566 | 79 | 66 | 7 | 1.237 | 0 | 0.0 |
| campylobacter | SP | high | 4 | 56 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| ecoli | DR | high | 4711 | 78111 | 23855 | 21990 | 376 | 0.481 | 66 | 0.084 |
| ecoli | SP | high | 860 | 12996 | 792 | 650 | 137 | 1.054 | 0 | 0.0 |
| efaecium | DR | high | 629 | 10281 | 2501 | 2352 | 47 | 0.457 | 35 | 0.34 |
| efaecium | SP | high | 21 | 203 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| enterobacter | DR | high | 422 | 6228 | 1828 | 1680 | 69 | 1.108 | 14 | 0.225 |
| enterobacter | SP | high | 51 | 465 | 23 | 19 | 3 | 0.645 | 0 | 0.0 |
| hpylori | DR | high | 92 | 1361 | 289 | 285 | 20 | 1.47 | 0 | 0.0 |
| hpylori | SP | high | 37 | 456 | 28 | 1 | 0 | 0.0 | 0 | 0.0 |
| kpneumoniae | DR | high | 1563 | 24831 | 7400 | 7047 | 57 | 0.23 | 23 | 0.093 |
| kpneumoniae | SP | high | 291 | 4405 | 234 | 192 | 43 | 0.976 | 0 | 0.0 |
| mtuberculosis | DR | high | 2840 | 64415 | 15450 | 13789 | 588 | 0.913 | 145 | 0.225 |
| mtuberculosis | SP | high | 410 | 234387 | 418 | 248 | 25 | 0.011 | 64 | 0.027 |
| ngonorrhoeae | DR | high | 46 | 825 | 212 | 149 | 6 | 0.727 | 4 | 0.485 |
| ngonorrhoeae | SP | high | 4 | 50 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |
| paeruginosa | DR | high | 2921 | 46265 | 13372 | 13090 | 173 | 0.374 | 12 | 0.026 |
| paeruginosa | SP | high | 597 | 9404 | 766 | 614 | 146 | 1.553 | 0 | 0.0 |
| pfalciparum | DR | high | 2710 | 75812 | 10741 | 10410 | 2859 | 3.771 | 2118 | 2.794 |
| pfalciparum | SP | high | 155 | 405323 | 239 | 71 | 176 | 0.043 | 0 | 0.0 |
| saureus | DR | high | 7709 | 135588 | 27882 | 25274 | 690 | 0.509 | 101 | 0.074 |
| saureus | SP | high | 1043 | 15722 | 1098 | 711 | 249 | 1.584 | 0 | 0.0 |
| smansoni | DR | high | 37 | 980 | 221 | 198 | 0 | 0.0 | 0 | 0.0 |
| smansoni | SP | high | 47 | 1936 | 15 | 8 | 6 | 0.31 | 0 | 0.0 |
| spneumoniae | DR | high | 1375 | 22589 | 4502 | 2693 | 141 | 0.624 | 90 | 0.398 |
| spneumoniae | SP | high | 92 | 832 | 0 | 0 | 0 | 0.0 | 0 | 0.0 |

## Verification: residual unproven actives in actual outputs

Actives whose qualifier does not prove a hit under the conservative rule (should be 0 post-fix; `n/a` = legacy output without `relation`).

| Pathogen | Actives | Unproven actives |
|---|---|---|
| abaumannii | 3909 | 0 |
| calbicans | 15678 | 0 |
| campylobacter | 193 | 0 |
| ecoli | 28345 | 0 |
| efaecium | 4590 | 0 |
| enterobacter | 3039 | 0 |
| hpylori | 733 | 0 |
| kpneumoniae | 8948 | 0 |
| mtuberculosis | 35625 | 0 |
| ngonorrhoeae | 487 | 0 |
| paeruginosa | 13707 | 0 |
| pfalciparum | 88147 | 0 |
| saureus | 63307 | 0 |
| smansoni | 993 | 0 |
| spneumoniae | 15299 | 0 |

## Worked examples (label flips across cutoff levels)

**abaumannii**

```
pair CHEMBL1953703 / MIC (DR, unit=umol.L-1, direction=-1)
molecule FYQMRZRXHKCOTG-VMTUJNLTSA-N: raw '>13.756175233038206' -> transformed value 55.69459866820291
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': False, 'high': False}
```

**abaumannii**

```
pair CHEMBL4051718 / MIC (DR, unit=umol.L-1, direction=-1)
molecule UNVQXEOCOPEOSI-NJAFHUGGSA-N: raw '>17.039476206501412' -> transformed value 24.14387334124045
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': False, 'high': False}
```

**calbicans**

```
pair CHEMBL2383960 / MIC (DR, unit=umol.L-1, direction=-1)
molecule VUKOWBDRNRSLFR-UJYXHFHISA-N: raw '>7.295648875010943' -> transformed value 11.12230756739562
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': True, 'high': False}
```

**calbicans**

```
pair CHEMBL889465 / MIC (DR, unit=umol.L-1, direction=-1)
molecule RNGIZOSWBIQUEI-MITPMGBZSA-N: raw '<12.470013514377149' -> transformed value 12.470013514377149
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': False, 'middle': True, 'high': True}
```

**campylobacter**

```
pair CHEMBL1120132 / MIC (DR, unit=umol.L-1, direction=-1)
molecule ZPTYJMSUNNGOPL-UHFFFAOYSA-N: raw '>13.8' -> transformed value 32.1
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': False, 'high': False}
```

**campylobacter**

```
pair CHEMBL1665006 / MIC (DR, unit=umol.L-1, direction=-1)
molecule OHKOGUYZJXTSFX-KZFFXBSXSA-N: raw '<10.404879888667786' -> transformed value 10.404879888667786
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': False, 'middle': True, 'high': True}
```

**ecoli**

```
pair CHEMBL4728694 / MIC (DR, unit=umol.L-1, direction=-1)
molecule WEPMVSQXDHVXOB-UHFFFAOYSA-N: raw '>16.05007623786213' -> transformed value 37.51826084101872
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': False, 'high': False}
```

**ecoli**

```
pair CHEMBL3366710 / MIC (DR, unit=umol.L-1, direction=-1)
molecule IQXLGSIFNVTMHP-ZJUUUORDSA-N: raw '>17.883329160556528' -> transformed value 195.5450917204253
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': False, 'high': False}
```

**efaecium**

```
pair CHEMBL3820038 / MIC (DR, unit=umol.L-1, direction=-1)
molecule NEZIIDJKGAEDTO-JLHYYAGUSA-N: raw '>19.58964589266343' -> transformed value 351.2546376588878
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': False, 'high': False}
```

**efaecium**

```
pair CHEMBL889463 / MIC (DR, unit=umol.L-1, direction=-1)
molecule IUICPAFYRDJMBJ-VTVXWXLCSA-N: raw '>12.48003194888179' -> transformed value 12.548680055323992
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': False, 'high': False}
```

**enterobacter**

```
pair CHEMBL978352 / IZ (SP, unit=mm, direction=1)
molecule CIDUJQMULVCIBT-MQDUPKMGSA-N: raw '>16.0' -> transformed value 16.0
cutoffs {'low': 10.0, 'middle': 20.0, 'high': 30.0}
ambiguous per level: {'low': False, 'middle': True, 'high': True}
```

**enterobacter**

```
pair CHEMBL1284229 / MIC (DR, unit=umol.L-1, direction=-1)
molecule NLVFBUXFDBBNBW-PBSUHMDJSA-N: raw '>17.111567419575636' -> transformed value 74.86080568942124
cutoffs {'low': 20.0, 'middle': 10.0, 'high': 5.0}
ambiguous per level: {'low': True, 'middle': False, 'high': False}
```


Per-dataset detail: `2026-07-01_censoring-interference.csv`