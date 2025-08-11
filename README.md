## Structure–Function Analyses for AD/FTD (MATLAB)

Code for analyses in “Functional network collapse in neurodegenerative disease” (Brown et al., 2025).

- Main scripts:
  - `adftd_structure_function_prep.m`: loads dataset, builds gradient codes (PCA on independent controls), reconstructs ROI-level FC, ComBat harmonization, vectorizes FC and gradient-pair covariances.
  - `adftd_structure_function.m`: runs manuscript figures (PLS: atrophy–FC; syndrome-specific FC contrasts; gradient–function links; coupled-oscillator/eigenmode dynamics).

### Data
- Large dataset: `adftd_structure_function_data.mat` (≈2.6 GB), placed in repo root.
- Hosted on Zenodo: [DOI: 10.5281/zenodo.16783268](https://doi.org/10.5281/zenodo.16783268).
- This repo is code-only; `*.mat` files are not tracked.

### Requirements
- MATLAB R2021b+ (tested on R2024a)
- Toolboxes: Statistics & Machine Learning (for `plsregress`, `classify`, `fitlm`, `regstats`)
- Optional: Mapping Toolbox (`wrapTo2Pi`); if missing, use `mod(x, 2*pi)`
- Included: `CircStat2012a/` (BSD) for circular stats
- R (optional, for cognition GAMs): `mgcv`, `gratia`, `ggplot2`, `dplyr`

### Quick start
1) Start MATLAB and set path:
```matlab
restoredefaultpath; rehash toolboxcache; addpath(genpath(pwd));
```
2) Optional preparation (variables are already in the `.mat` for reproducibility):
```matlab
adftd_structure_function_prep
```
3) Run analyses and figures:
```matlab
adftd_structure_function
```

### What the scripts do
- `adftd_structure_function_prep.m`
  - PCA on independent controls → gradient codes; project all scans into this basis.
  - Compute gradient covariance; reconstruct ROI-level FC via gradient weights.
  - ComBat harmonization across scanners for atrophy and gradient covariance.
  - Vectorize FC upper triangle (246 ROIs → 30,135 edges) and 6×6 gradient-pair covariances.

- `adftd_structure_function.m`
  - PLS (Figure 1): atrophy (321×246) vs FC (321×30,135); sign aligned so higher scores = higher atrophy; plots score correlations and FC loadings.
  - Syndrome FC (Figure 2): LDA selects typical cases; observed vs reconstructed FC contrasts.
  - Gradient–function link (Figure 4): regress PLS `Yscores` on gradient-pair covariances; report R², significant pairs, contributions.
  - Dynamics / eigenmodes (Figures 5–6; Supp): estimate coupling from gradient 1st/2nd derivatives (`latent_derivatives`), simulate with `coupled_oscillator_function`, compare predicted vs observed FC, eigenmode magnitudes/angles vs gradient covariance.

### Optional: Cognition GAMs (R)
- Script: `cogniton_gam.R` fits GAMs for multiple cognitive measures using structural and functional PLS scores with covariates (age, sex, education, scanner, motion).
- Inputs live in `cognition_data/` as `yX_*.txt`. Update the hard-coded `data_dir` in the script if needed to point to the repo path.

### Key files
- Main: `adftd_structure_function_prep.m`, `adftd_structure_function.m`
- Harmonization: `combat.m`, `combat_scripts/`
- Dynamics: `latent_derivatives.m`, `coupling_parameters.m`, `coupled_oscillator_function.m`
- Plotting: `plot_matrix.m`, `plotcorr.m`
- Ancillary: `CircStat2012a/` (circular stats), `gradient_maps/` (NIfTI maps for display), `cognition_data/` (inputs for GAMs), `cogniton_gam.R`

### Citation
If you use this code:
- Paper: Brown J, et al. Functional network collapse in neurodegenerative disease. 2025.
- Code: https://github.com/jbrown81/structure_function
- Data: https://doi.org/10.5281/zenodo.16783268

### License & patent
Apache License 2.0.

### Contact
Jesse Brown — jesse@radiata.ai or jbrown81@gmail.com