# ADNI Cortical Thickness JSD

Code for the ADNI analyses in *Longitudinal Changes in Cortical Thickness Distribution Networks Across the Clinical Spectrum of Alzheimer's Disease: Associations With Cognition*.

The workflow starts with quality-controlled FreeSurfer reconstruction outputs and ends with the ADNI statistical results, Figure 1–3 and the ADNI portions of Tables 1–2. It does not run MRI conversion, cortical reconstruction, FastSurfer, FreeSurfer reconstruction commands or Docker preprocessing. It does not reproduce the external cohorts.

The repository contains code, numerical configuration, ROI/edge dictionaries and aggregate validation records. ADNI images, reconstructed surfaces, clinical records, individual JSD features, fold assignments and trained model parameters are not distributed. Obtain study data through authorized [ADNI access](https://adni.loni.usc.edu/data-samples/adni-data/).

## Implementation

The JSD numerical functions and configuration come from the supplied `ADNI_JSD256_h01_range6_transfer_package.zip`. The statistical implementations come from the formal ADNI analysis that produced the paper. New command-line entry points connect these components and replace computer-specific paths. The source mapping is in `docs/version_audit.csv`.

The numerical extraction is area weighted, uses reflected Gaussian KDE with bandwidth 0.1 mm, 256 equal-width bins over 0–6 mm and a separate upper-tail bin. The 257 probability masses retain the supplied integration, h/32 linear binning, sparse-center fallback, double precision, normalization and base-2 JSD implementation. ROI order is fixed at 62 parcels; edge order is the R column-major upper triangle of their 62-by-62 matrix. No clipping of the upper tail or replacement of the numerical method is introduced.

The classifier uses the original two-block ridge objective, parameter grid, inner log-loss selection and numerical solver fallback. Standardization is calculated separately in each inner and outer training split using the population SD. The outer design has five repeats of five folds and four inner folds. Generated patient-level fold files remain local. Full-classifier AUC is evaluated from repeat-averaged out-of-fold probabilities. Its conditional 95% CI follows the recorded class-stratified Python bootstrap with 2,000 draws. The reported AUC of 0.925 is the full M3 classifier's performance, not that of the JSD contribution score alone.

Fixed scoring extracts only the JSD contribution. Each source model uses its training-CN mean and population SD. A classifier participant receives the held-out fold in each repeat; other participants receive the prescribed hash-based fold. The same source models are used throughout each participant's visits. The five contributions are averaged and normalized with the mean and sample SD of all retained ADNI baseline scores. The program derives these constants from the inputs; it does not substitute published values.

Longitudinal models use the original REML mixed models, diagonal random intercept/time-slope covariance and recorded intercept-only/optimizer fallback. Inference uses participant-clustered CR1 covariance and participant N minus one degrees of freedom. Baseline cognitive regressions use HC3 covariance. Between/within components are calculated in each fitted analysis sample. Cognitive observations are MRI-linked records retained by the formal complete-case and distinct-date rules; no additional assessment-date deduplication is applied.

Connection models test every one of the 1,891 edges. BH correction is performed separately by model and contrast over the full 1,891-edge family. The primary discovery set is L0 AD–CN q < 0.05. Group slope signs determine the convergence/divergence categories. Homotopic results retain whole-network correction. Score slope differences have a separate three-comparison BH family. Nonsignificant effects and confidence intervals crossing zero remain in the output.

## Environment

Two environment files reflect the installed runtimes used for validation:

- `environment.yml`: Python 3.12, NumPy 2.3.5, pandas 3.0.1 and openpyxl 3.1.5.
- `environment-figures.yml`: Python 3.10, NumPy 1.24.3, pandas 2.3.3, SciPy 1.15.3, Matplotlib 3.7.1, nibabel 5.2.0 and Pillow 12.1.1.

R is installed separately. The validated R runtime was R 4.6.0 with glmnet 5.1, nlme 3.1-169 and Matrix 1.7-5; versions are also recorded in `docs/software_versions.csv`. Install these R packages into an accessible R library and set `R_LIBS_USER` if it is not a default library. The Python environment files specify dependencies; a fresh Conda solve was not performed during packaging. The validation report identifies the actual interpreters used. Check an installed environment from the repository root:

To create the Python environments with Conda, use `conda env create -f environment.yml` and `conda env create -f environment-figures.yml`. Activate `adni-jsd-analysis` for analysis commands and `adni-jsd-figures` for plotting. Ensure `Rscript` is available on the terminal PATH.

```text
python scripts/check_environment.py
Rscript --vanilla scripts/check_environment.R
python scripts/check_environment.py --figures
```

Use the figure environment for plotting commands. Arial was used for the manuscript figures. Another font installation can change text metrics and image appearance.

## Local inputs and paths

Run commands from this repository's root. Keep inputs and generated analysis files outside Git tracking. The defaults are `data/` and `analysis/`, both ignored by `.gitignore`. Set `ADNI_DATA_DIR` and `ADNI_WORK_DIR` to use other authorized local locations. Set `ADNI_FS_ROOT` to the directory containing reconstructed scan folders.

For example, after placing the authorized inputs in sibling directories, set the paths in PowerShell:

```powershell
$env:ADNI_DATA_DIR = (Resolve-Path "../private_inputs").Path
$env:ADNI_FS_ROOT = (Resolve-Path "../freesurfer_subjects").Path
$env:ADNI_WORK_DIR = Join-Path (Get-Location) "analysis"
$env:ADNI_ANATOMY_DIR = (Resolve-Path "../display_surface").Path
```

Replace the relative input-directory names with the actual local locations. The display reference is needed for Figure 3; it is independent of the reconstruction folders used to calculate JSD.

Required input tables are described by the header-only files in `docs/schemas/`:

- `scan_manifest.csv`: one row per reviewed scan, with `PTID`, `MRI_IMAGE_ID`, `MRI_DATE` and a relative `FOLDER`. Row order determines the extraction cache index. The manifest must be supplied locally; no participant roster is included.
- `clinical_mri.csv`: the study's already curated, MRI-centered clinical table with unique `PTID`/`MRI_IMAGE_ID` keys. It supplies baseline diagnosis, MRI age, sex, education, ADAS13, MMSE and the clinical assessment date used in the formal models. The matcher retains the exact MRI keys and the original eligibility rule of at least two valid MRI records per participant. It does not rebuild the underlying ADNI clinical database or impose a new nearest-date/window rule.

Model-building fields are `DX.bl`, `AGE_AT_MRI_EST`, `PTGENDER`, `PTEDUCAT`, `ADAS13`, `MMSE` and `EXAMDATE`. The schema also retains the ancillary metadata columns required by the formal input builder; ancillary values may be missing. `CORTEX_CM3` and global mean cortical thickness come from the reconstruction outputs. All ROI means, ROI SDs and JSD features come from the extraction steps.

Each reconstructed scan directory must retain:

```text
surf/lh.white                         surf/rh.white
surf/lh.thickness                     surf/rh.thickness
label/lh.aparc.DKTatlas.annot          label/rh.aparc.DKTatlas.annot
label/lh.cortex.label                 label/rh.cortex.label
stats/aseg.stats
stats/lh.aparc.DKTatlas.mapped.stats   stats/rh.aparc.DKTatlas.mapped.stats
```

The original DKT `.stats` filenames are accepted when `.mapped.stats` is absent. Retain the existing `sulc`, `curv`, `pial` and reconstruction logs as well. Available sulc/curvature measurements and the logged white-surface curvature check participate in the original common-valid mask. Those logs are read locally to preserve that numerical rule; they are not distributed or collected as a public MRI processing history. ROI QC requires at least 30 retained vertices, effective area-weighted N of at least 20 and at least 90% retained cortical area.

Figure 3 additionally requires an authorized surface reference directory through `ADNI_ANATOMY_DIR`, containing `surf/{lh,rh}.pial` and `label/{lh,rh}.aparc.DKTatlas.annot`. The manuscript used a study-specific display surface, which is not included. Using another reference changes the surface rendering but not the connection statistics or ROI endpoint counts.

## Run order

The commands below run the study workflow; they require the authorized inputs described above. This delivery did not execute all scans or retrain all 25 classifier checkpoints.

1. Extract the reviewed ROI vertices and white-surface area weights:

   ```text
   Rscript --vanilla scripts/01_roi_extraction/extract_rois.R 1 1
   ```

   Output: `roi/scan_*.rds`, including the local raw ROI arrays, QC and cortical statistics. For parallel extraction, run distinct worker IDs with the same worker count; row assignment is fixed.

2. Calculate probability masses and JSD features:

   ```text
   Rscript --vanilla scripts/02_jsd_calculation/calculate_jsd.R
   ```

   Output: `jsd/JSD_features.csv`, `ROI_QC.csv`, `cortical_stats.csv`, `edge_order.csv`, and per-scan RDS files containing the 257-bin masses and 62-by-62 matrices. A failed ROI/scan stops completion; review the input manifest rather than silently replacing or imputing data.

3. Match clinical records and construct the formal analysis inputs and folds:

   ```text
   python scripts/03_data_preparation/match_clinical.py
   python scripts/03_data_preparation/prepare_model_inputs.py
   ```

   Output: `data/clinical_jsd.csv`; `score/inputs/all_scans.csv`, `classifier_input.csv`; local outer/inner fold files and numerical protocol. The input builder also accepts an explicit path to a matching curated master CSV as its first argument.

4. Fit the source classifier and derive/apply the fixed score:

   ```text
   Rscript --vanilla scripts/04_score_derivation/run_nested_cv.R
   Rscript --vanilla scripts/04_score_derivation/export_parameters.R
   python scripts/04_score_derivation/apply_score.py
   python scripts/04_score_derivation/evaluate_oof.py
   ```

   Output: local classifier checkpoints, OOF predictions, selected parameters, `M3_parameters.csv`, `scored_MRI.csv`, `baseline_scores.csv`, `score_scale.json`, `classifier_performance.csv` and `repeat_AUC.csv`. These participant-level and fitted-model outputs remain local. Verified same-protocol checkpoints are reusable; incompatible checkpoints are rejected. The original worker interface is also retained: `worker WORKER_ID WORKER_COUNT bounded`, followed by a `collect` invocation after all workers finish.

5. Fit longitudinal score models:

   ```text
   Rscript --vanilla scripts/05_longitudinal_models/fit_trajectories.R
   ```

   Output: `core/results/trajectory_results.csv`, coefficient tables, fitted model objects and local inclusion records for L0 and global-MRI-adjusted L1.

6. Fit cognitive associations:

   ```text
   Rscript --vanilla scripts/06_cognitive_models/fit_longitudinal_cognition.R
   Rscript --vanilla scripts/06_cognitive_models/fit_baseline.R
   ```

   Output: longitudinal ADAS13/MMSE between/within estimates; C0/C1/C2 baseline regression coefficients and the prespecified MCI-only C1 coefficients. Sample selection and adjustment follow the formal scripts.

7. Fit and collect all connection models:

   ```text
   python scripts/07_connection_models/prepare_edges.py
   Rscript --vanilla scripts/07_connection_models/fit_edges.R 1 1
   Rscript --vanilla scripts/07_connection_models/collect_edges.R
   ```

   Output: local column-major float64 features, model checkpoints, `edges/results/all_edge_effects.csv`, separate L0/L1 tables, fixed coefficients and whole-network FDR results. All valid edge/model fits are required by the collector. Worker IDs/counts can partition the fixed edge order without changing the model.

8. Export the figure inputs, predictions and final ADNI tables/figures:

   ```text
   python scripts/08_figures_tables/prepare_figures.py
   Rscript --vanilla scripts/08_figures_tables/export_curves.R
   python scripts/08_figures_tables/render_figures.py
   python scripts/08_figures_tables/tables.py
   ```

   Output: `figures/figure1` through `figure3`, each with vector PDF/SVG and 600 dpi PNG/TIFF; `tables/Table1_ADNI.csv`, `Table2_ADNI.csv` and a printable HTML table. Figure 2 uses participant means or within-person deviations for its descriptive scatters and saved adjusted mixed-model coefficients for its annotations. Figure 3B masks by primary whole-network FDR; the complete unmasked edge results remain in CSV. Edge plots display 10^-3 JSD/year; CSV estimates retain raw JSD/year.

## Validation and publication checks

```text
Rscript --vanilla tests/test_jsd.R
python tests/verify_study_results.py
python tests/check_public_release.py
```

`test_jsd.R` checks the numerical configuration, ROI/edge order, mesh area calculation, probability normalization, upper tail and dense/sparse agreement using mathematical test inputs. These inputs are not simulated study results and are not used for figures. `verify_study_results.py` checks generated outputs against published aggregate reference values in `config/expected_results.json`; the reference values never enter model fitting or sample selection. `check_public_release.py` scans the candidate source tree for prohibited artifacts and detectable identifiers, personal paths and credentials.

See `docs/validation_report.md` for the tests actually executed, `docs/test_results.csv` for their aggregate results, and `docs/reproduction_limits.md` for remaining limits. Do not infer full cohort-level extraction or full classifier retraining from a single-scan numerical check or a checkpoint audit.

Before distributing generated outputs, review their content under the applicable ADNI agreement. `.gitignore` is a convenience, not a data-use authorization. The software licensing status is recorded in `LICENSE_NOTICE.md`. This repository has not been uploaded automatically.
