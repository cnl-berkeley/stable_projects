# Overlooked neuroanatomical correlates of face processing and developmental prosopagnosia in posteromedial cortex

Data and code for the manuscript *Overlooked neuroanatomical correlates of face processing and developmental
prosopagnosia in posteromedial cortex*.

For questions and/or additional data requests please email Joe Kelly (josephkelly1@northwestern.edu) and/or
Kevin Weiner (kweiner@berkeley.edu).

## Repository structure

```
data/                         processed data used by the analyses (data dictionary below)
PMC_DP_analysis_vF.R          face selectivity and sulcal incidence           Figures 2, 3; Table 1
demographics_covariates.R     demographics, group matching, covariates, CFMT  Methods; Figure S1
PMC_face_DP_vF.ipynb          LASSO analysis of sulcal depth and CFMT          Figure 6; Table S1
power/lasso_power_simulation.py   power analysis for the LASSO analysis        Methods, "Power analysis"
centroid/centroid_distance_stats.R  face-patch to sulcus distances             Figure 4
dice/dice_prcusp_vs_other_sulci.py  sulcus overlap with the face-selective region (within participants)  Figure 5
dice/dice_group_tests.R             the same overlap, NT vs DP                  Figure 5
surface_producers/            scripts that computed the surface-based data files from FreeSurfer
                              surfaces (documentation only)
outputs/                      outputs of every script, as used in the manuscript
```

## What each script reproduces

| Script | Manuscript content |
|---|---|
| `PMC_DP_analysis_vF.R` | Results: sspls-v incidence and face selectivity (all statistics); sulcal–functional gradient (all statistics); HCP replication. Supplemental Material: face selectivity of all PMC sulci. Figure 2B, 2C; Figure 3B; Table 1; total number of labelled sulci |
| `demographics_covariates.R` | Abstract (percent female); Methods: participants of Datasets 1 and 2, group matching (age, sex), CFMT descriptives, LASSO sample, logistic regressions of sulcal presence on age and sex; participants of Dataset 3 (HCP; computed only if `data/hcp_demographics.csv` is supplied, see below); Figure S1 |
| `PMC_face_DP_vF.ipynb` | Results: right-hemisphere LASSO (α, RMSE<sub>CV</sub>, coefficients, fold selection), bootstrap summary, AIC comparison, cross-validated adjusted R², control models (left hemisphere, cortical thickness, DPs); Figure 6A, 6B; Table S1 |
| `power/lasso_power_simulation.py` | Methods, "Power analysis" (effect detectable at 80% selection, multicollinearity, selection and R² at β = .30 and .50, null selection rate) |
| `centroid/centroid_distance_stats.R` | Methods: fragmentation of the face-selective region. Results: distances to the precuneal sulci, gradient model, group effects at 2.5/5/10%, prculs-d positional control (position, distance, within-group, between-group). Figure 4 |
| `dice/dice_prcusp_vs_other_sulci.py` | Results and Supplemental Material: prcus-p vs every other consistently present sulcus at 0, 0–5 and 0–10 mm, at 2.5/5/10%; highest mean DSC at each distance; DSC peak; mcgs non-overlap; area of the face-selective region; Figure 5 caption group means |
| `dice/dice_group_tests.R` | Results and Supplemental Material: NT vs DP permutation tests (group, sulcus-by-group) and prcus-p rank-sum tests, uncorrected and Holm-corrected across the six window-by-hemisphere tests, at 2.5/5/10%; Figure 5B |

Figure 1, Figure 2A, Figure 3A, Figure 5A, Figure 6C and Figures S2–S6 are anatomical renderings or literature
panels and are not produced by these scripts.

## Requirements

* **R** 4.2.3 with tidyverse 2.0.0 (dplyr, tidyr, ggplot2), rstatix 0.7.2, lme4 1.1-32, lmerTest 3.1-3,
  emmeans 1.8.5, nlme 3.1-162, patchwork 1.3.0.
* **Python** 3.10 with numpy 2.2, pandas 2.3, scipy 1.15, scikit-learn 1.7.2, Jupyter (nbconvert 7, ipykernel 7)
  and rpy2 3.5 (for the notebook's figures, using an R installation with ggplot2 and dplyr).
  The notebook works with scikit-learn versions before and after the removal of `mean_squared_error(squared=False)`.
  `power/lasso_power_simulation.py` calls scikit-learn's internal coordinate-descent routine for speed and checks it
  against the public `lasso_path` before running; it was tested with scikit-learn 1.7.2.
* The scripts in `surface_producers/` additionally need nibabel, per-participant FreeSurfer surfaces and labels
  (not shared), an identifier map from FreeSurfer subject names to the participant labels used here (not shared;
  see `surface_producers/surface_io.py`), and, for the sulcal centroids, the CNL scalpel package
  (https://github.com/b-parker/CNL_scalpel).

## How to run

Run every script from the repository root. The scripts are independent of one another (each reads only `data/`)
and can be run in any order.

```
Rscript PMC_DP_analysis_vF.R                      # seconds
Rscript demographics_covariates.R                 # seconds
Rscript centroid/centroid_distance_stats.R        # seconds
python3 dice/dice_prcusp_vs_other_sulci.py        # seconds
Rscript dice/dice_group_tests.R                   # under a minute (10,000 permutations per test)
jupyter nbconvert --to notebook --execute PMC_face_DP_vF.ipynb   # about 5 minutes
python3 power/lasso_power_simulation.py           # 1-2 hours on 8 cores
```

Options for the resampling steps (defaults equal the values used in the paper):

* `Rscript dice/dice_group_tests.R N_PERM` — number of permutations (default 10,000). The random seed is fixed,
  so the default reproduces the reported p-values exactly.
* `python3 power/lasso_power_simulation.py --iterations N --betas ...` — simulated datasets per cell (default 150)
  and planted effects (default .10–.70). Every simulated dataset has a fixed seed, so the default reproduces the saved
  full run (`outputs/power/saved_full_run/`). `--from-saved` summarizes the saved full run without simulating.
* LASSO bootstrap: set `RUN_BOOTSTRAP=1` (and optionally `N_BOOTSTRAP`, default 2,000) before running the notebook.
  The bootstrap in the paper (Table S1) was run once without a fixed random seed; its saved results are in
  `outputs/lasso/bootstrap/`, and a re-run will give slightly different values. The original run took about
  9.5 hours. Re-run results are written to `outputs/lasso/bootstrap/rerun/` and never overwrite the saved results.

## Data dictionary

All files are comma-separated, one header row. `sub` is the participant identifier: for Datasets 1 and 2 the
label used in Figures S2–S5 (`NT1`–`NT43`, `DP1`–`DP39`), for HCP participants the public HCP subject ID; `group` is `Controls`
(neurotypical, NT), `DPs` (developmental prosopagnosia) or `HCP_Controls` (Human Connectome Project); `hemi` is `lh`
or `rh`. Sulcus abbreviations follow the manuscript, except that the splenial sulcus (spls) is labelled `sbps`
(its alternative name, the subparietal sulcus). Datasets 1 and 2 are pooled for the structural analyses
(NT n = 43, DP n = 39); Dataset 1 alone (NT n = 25, DP n = 22) has functional data.

| File | Rows | Columns |
|---|---|---|
| `demographics.csv` | one per participant of Datasets 1 and 2 (82) | `sub`; `dataset` (1 or 2); `group`; `age` (years); `sex` (F/M); `CFMT` (Cambridge Face Memory Test score out of 72; missing for one NT) |
| `sulcus_presence.csv` | participant × hemisphere × sulcus for the 12 PMC sulci (1,968) | `sub`, `group`, `hemi`, `sulcus`, `present` (1 = the sulcus was identified in that hemisphere) |
| `face_selectivity.csv` | participant × hemisphere × sulcus, sulci present only | `sub`, `group`, `hemi`, `sulcus`, `mean_face_activation`: mean face selectivity across the sulcus; Dataset 1: faces − objects, % signal change (leave-one-run-out, Jiahui et al., 2018); HCP: faces − all other categories, z (Chen et al., 2023) |
| `sulcal_depth.csv` | participant × hemisphere, Datasets 1 and 2 (162; one NT without a CFMT score has no row) | `sub`, `group`, `hemi`, `CFMT`, and the normalized mean sulcal depth (FreeSurfer .sulc divided by the hemisphere's maximum) of the nine sulci in the LASSO analysis: `mcgs`, `pos`, `sbps`, `prculs-d`, `prcus-p`, `prcus-i`, `prcus-a`, `ifrms`, `sspls-v` (empty when the sulcus is absent) |
| `sulcal_thickness.csv` | as `sulcal_depth.csv` | as `sulcal_depth.csv`, with mean cortical thickness (mm; FreeSurfer mris_anatomical_stats) |
| `centroid_distance.csv` | participant × hemisphere × threshold × sulcus, Dataset 1 (1,692) | `threshold_pct` (face-selective region = top 2.5, 5 or 10% of PMC vertices); `sulcus` (prcus-a/i/p, prculs-d, and mcgs and pos, shown for illustration in Figure 4D); `geodesic_mm`: geodesic distance on the native white surface from the area-weighted centroid of the face-selective region to the area-weighted centroid of the sulcus |
| `face_patch_geometry.csv` | participant × hemisphere × threshold, Dataset 1 (282) | `n_pmc_vertices` (vertices in the PMC territory: six Destrieux parcels); `n_suprathreshold_vertices`; `n_components` (connected components of the suprathreshold vertices); `top_component_mass_fraction` (share of total mass, area × t, in the highest-mass component); `top_component_vertices`; `top_component_area_cm2` (area of the face-selective region used for the centroid and Dice analyses) |
| `sulcal_ap_position.csv` | participant × hemisphere × sulcus (prcus-p, prculs-d), Dataset 1 (188) | `ap_position`: mean anterior–posterior position of the sulcal label, 0 = pos, 1 = mcgs (geodesic) |
| `dice_overlap.csv` | participant × hemisphere × threshold × sulcus × distance, Dataset 1 (11,280) | `threshold_pct`; `sulcus` (the eight sulci present in every hemisphere); `radius_mm` (dilation of the sulcal label, 0–10 mm in 2.5-mm steps); `dice`: Dice similarity coefficient between the dilated label and the face-selective region (surface area, native surface) |

Provenance of the surface-based files: `centroid_distance.csv` from `surface_producers/centroid_distance.py`,
`face_patch_geometry.csv` from `surface_producers/patch_geometry.py` (with the region area from
`dice_dilated_native.py`), `sulcal_ap_position.csv` from `surface_producers/ap_position.py`, and
`dice_overlap.csv` from `surface_producers/dice_dilated_native.py` (run at 5, 2.5 and 10%). Sulcal labels,
depth and thickness were extracted from manually defined labels on each participant's FreeSurfer surfaces
(see Methods).

**HCP participants.** Face selectivity for the 71 HCP participants is in `face_selectivity.csv`. Their age and sex are not
redistributed here (exact age is HCP Restricted Data); to reproduce the Dataset 3 demographics, place them in
`data/hcp_demographics.csv` (columns `sub`, HCP subject ID as in `face_selectivity.csv`; `age`, years; `sex`, F/M;
one row per participant) and re-run `demographics_covariates.R`. HCP data are from the WU-Minn HCP Young Adult S1200 release
(https://www.humanconnectome.org).


## Notes on the statistics

* Omnibus ANOVAs use `rstatix::anova_test` as described in the Methods; post hoc tests use a linear mixed-effects
  model fit to every available observation. One-sample tests of face selectivity against zero are Bonferroni-adjusted
  across the two hemispheres within group for the sspls-v, and across the 12 sulci within each hemisphere and group
  for the tests of all PMC sulci; `PMC_DP_analysis_vF.R` writes both families.
* The AIC comparison and the cross-validated R² in the notebook use the held-out predictions of the fully nested
  LASSO at the optimal α (1.35) for the selected model, and leave-one-out ordinary least squares for the full model.
