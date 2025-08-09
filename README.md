# Structure–Function Analyses for AD/FTD (MATLAB)
Code for the analyses in *Functional network collapse in neurodegenerative disease* (Brown et al., 2025).

This repo has two main scripts:

- **`adftd_structure_function_prep.m`** — loads the dataset, builds gradient codes (PCA), runs ComBat harmonization, reconstructs ROI-level FC, and prepares all variables.
- **`adftd_structure_function.m`** — runs the manuscript analyses/figures (PLS on atrophy vs FC, syndrome-specific FC contrasts, gradient–function links, coupled-oscillator/eigenmode analyses).

> **Data policy:** This repo is **code-only**. The large dataset (`adftd_structure_function_data.mat`, ~2.6 GB) is hosted on Zenodo.  
> **Zenodo DOI:** https://doi.org/10.5281/zenodo.16783268. Place the `.mat` next to the `.m` files.

---

## Requirements
- **MATLAB** R2021b or newer (tested on R2024a)
- **Toolboxes:** Statistics & Machine Learning (for `plsregress`, `classify`, `fitlm`, `regstats`)
- **Optional:** Mapping Toolbox (for `wrapTo2Pi`) — if missing, use `mod(x, 2*pi)` instead
- **Included:** `CircStat2012a/` (BSD) for circular stats utilities

Runs on macOS/Windows/Linux. The full dataset needs a few GB of free RAM/disk.

---

## Get the code
```bash
git clone https://github.com/jbrown81/structure_function.git
cd structure_function
```

---

## Get the data
- **Full dataset (2.6 GB):** `adftd_structure_function_data.mat` → **Zenodo DOI: https://doi.org/10.5281/zenodo.16783268**  
  Put the file in the repo root.
---

## Quick start (full pipeline)
1) **Start MATLAB**, clear path, and add subfolders to the path:
```matlab
restoredefaultpath; rehash toolboxcache; addpath(genpath(pwd));
```

2) **Optional: Prepare inputs** (loads the `.mat`, builds gradients, ComBat, FC vectors; running this script is not necessary for running the main analysis script as all the variables created by this script are also already preloaded in adftd_structure_function_data.mat):
```matlab
adftd_structure_function_prep
```

3) **Run analyses & figures**:
```matlab
adftd_structure_function
```

Figures display in MATLAB as they’re generated.

---

## What each script does

### `adftd_structure_function_prep.m`
- Loads `adftd_structure_function_data.mat`
- PCA on independent controls → gradient codes; projects all subjects
- Computes gradient covariance; reconstructs ROI-level FC from gradient covariance + ROI weights
- ComBat harmonization (scanner batch) for atrophy & gradient covariance
- Vectorizes FC upper triangle (`pls_fc`) and gradient-pair covariances (`pls_grad_cov_n`)
- Leaves in memory (examples):  
  `atrophy_300_harmonized`, `pls_fc`, `pls_fc_pairs`, `keep_pls_fc`,  
  `grad_covs_all_harmonized`, `pls_grad_cov_n`, `pls_grad_pairs_n`,  
  `fc_mats_via_codes`, `roi_comp_slopes_fc_in_grp2`, `codes_fc_in_all_proj_grp2`,  
  `dxs_all`, `keep_str_fc`, `scanners_all`, etc.

> See the header of `adftd_structure_function_prep.m` for the full list of **required external variables** expected inside the dataset.

### `adftd_structure_function.m`
- **PLS (Figure 1):** atrophy (321×246) vs FC (321×30135); flips signs so higher scores = higher atrophy; plots score correlations and FC loading matrices
- **Syndrome-specific FC (Figure 2):** LDA to select typical cases; observed vs reconstructed FC contrasts
    
- **Gradient–function link (Figure 4):** regress `Yscores` on gradient pair covariances; report R², significant pairs, term contributions
- **Dynamics / eigenmodes (Figures 5–6; Supp Fig 10):** estimate coupling; simulate; compare predicted vs real FC; compute eigenmode magnitudes/angles and correlate with gradient covariance

---

## Repo contents (key files)
- `adftd_structure_function_prep.m`, `adftd_structure_function.m` — main scripts  
- `combat.m`, `combat_scripts/` — ComBat harmonization  
- `coupling_parameters.m`, `latent_derivatives.m`, `coupled_oscillator_function.m` — dynamical systems helpers  
- `plot_matrix.m`, `plotcorr.m` — plotting helpers  
- `CircStat2012a/` — circular stats utilities (BSD)  
- `gradient_maps/`, `cognition_data/`, `cogniton_gam.R` — ancillary resources

> **Not tracked:** `*.mat` files (large data). This repo intentionally avoids Git LFS.

---

## Troubleshooting
- **`wrapTo2Pi` not found:** Replace calls with `mod(x, 2*pi)` or install Mapping Toolbox.  
- **Out-of-memory with full dataset:** Try the demo dataset or run sections individually.  
- **Path issues:** Re-run `restoredefaultpath; rehash toolboxcache; addpath(genpath(pwd));`.

---

## Citation
If you use this code, please cite:
```
Brown J, et al. Functional network collapse in neurodegenerative disease. 2025.
Code: https://github.com/jbrown81/structure_function
```

---

## License & patent
Licensed under **Apache License 2.0** (see `LICENSE`).  
This implementation practices Radiata’s patented methods (U.S. Patent **[fill in number]**).  
Radiata grants Apache-2.0 licensees a **royalty-free, worldwide** patent license **solely** for use with software released under Apache-2.0 or a compatible OSI license.

---

## Contact
Jesse Brown — `jesse@radiata.ai`
