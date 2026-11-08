# Validation report

This report describes execution of the packaged code on authorized local inputs. Participant-level inputs, checkpoints, study surfaces and detailed working logs are outside the release.

## Source checks

- JSD authority: the supplied `ADNI_JSD256_h01_range6_transfer_package.zip` (SHA-256 `81afd7c88f1ba0fc9390e3f8e875da4ca79d9f3e43666072603332af0e6993f7`). Its numerical configuration was preserved. Thirty-five used function signatures/bodies match the supplied package.
- Sixteen formal statistical function expressions are unchanged, including the ridge fitting/preprocessing/selection functions, REML/CR1 functions, HC3 function and connection fitting function. Public changes concern paths, interfaces, output names, observed-count reporting and an input file-length check.
- The formal model-input builder was found and integrated. The clinical matching entry point uses the formal exact-key join and valid-MRI participant count rule. The MRI-centered clinical table is an authorized input; rebuilding the underlying clinical database is outside this repository.
- All required numerical and formal statistical source components within the agreed scope were found. New entry points are explicitly distinguished from the source computation in `version_audit.csv`.

## Executed checks

| Component | Actual execution |
|---|---|
| FreeSurfer to ROI to JSD | One real retained MRI. Both new workflow stages were executed. Its 1,891 JSD features agree with the formal input within approximately 5e-16; ROI moments and global MRI measures also agree within floating-point tolerance. |
| JSD numerical contracts | Fixed configuration/order, triangle area allocation, base-2 identities, reflected mass normalization, nonzero upper-tail retention and dense/sparse equivalence passed. Mathematical test inputs were not used as study data. |
| Clinical matching/input construction | Existing valid JSD features were joined to the authorized clinical table, yielding 770 people/3,593 MRI. The 352-row classifier input agrees numerically with the formal input. Model-relevant all-scan values agree within 5.33e-15. |
| Cross-validation | All 1,760 outer and 7,040 inner fold records regenerated exactly. All 25 saved checkpoints were inspected, including 100 inner partitions. Train/validation/test isolation, training-only population-SD standardization and six-model held-out probabilities passed. No complete repeated-nested-CV retraining was performed. |
| Fixed scoring | The parameter exporter was executed on all 25 formal checkpoints; all 47,425 exported coefficient rows match. Fixed source parameters were applied to all 3,593 existing JSD records. The assigned model rules, training-CN reference normalization and final baseline normalization reproduce the formal score within 1e-10. |
| Classifier evaluation | All six AUC/95% CI triples from saved OOF predictions and all 30 repeat AUC values agree. M3 is 0.925 (0.898–0.949). The recorded conditional bootstrap was executed without refitting classifiers. |
| Score/cognitive models | Both trajectory models, both longitudinal cognitive models and eight baseline/subgroup regressions were executed with the formal specifications. Maximum estimate/interval discrepancy was 3.22e-11. The primary analysis counts agree. |
| Connection models | One connection was refitted under L0 and L1; coefficients, CR1 covariance, contrasts and CIs match exactly. All 1,891 saved formal connection checkpoints were collected with the packaged collector. All 11,346 group-slope/contrast rows and their FDR results agree with the formal output; 277 primary discoveries and four homotopic discoveries were retained. The complete 3,782-fit batch was not refitted. |
| Figures/tables | Figure 1–3 were generated with the formal local surface reference. Their 600 dpi PNG pixels match the formal figures exactly; PDF renders were reviewed and no visible text falls outside the canvas. All nine reported ADNI Table 1 entries and eight ADNI Table 2 rows agree. |

Floating-point tolerances distinguish numerical reproduction from equality of formatted paper values. The largest model discrepancy is many orders below the manuscript display precision. PDF/SVG byte identity is not claimed because export metadata and object identifiers can differ.

The all-scan input audit checks the fields actually used by the formal models. The ancillary visit-level `DX` column in the public clinical join is retained from the supplied clinical table; it is not the grouping variable. Formal grouping remains derived from `DX.bl` exactly as in the source analysis.

## Environment and limits

Actual execution used Windows x86-64, R 4.6.0, the installed analysis Python 3.12.14 and figure Python 3.10.6 runtimes. Package versions are in `software_versions.csv`. The installed glmnet 5.1 package reports that it was built under R 4.6.1; it loaded and its saved-model calculations were checked under the actual R 4.6.0 runtime. A fresh Conda installation/solve and other operating systems were not tested.

This is a tested source-code integration, not a claim that every MRI or every trained model was independently recomputed. See `test_results.csv` and `reproduction_limits.md`. No history of individual MRI processing runs is distributed.
