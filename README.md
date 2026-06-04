# Efficient Cross-Fitting Estimation of Semiparametric Spatial Point Processes

Replication code for the paper:

> **Efficient, Cross-Fitting Estimation of Semiparametric Spatial Point Processes**
> [[arXiv:2410.04359]](https://arxiv.org/pdf/2410.04359)

---

## Overview

This repository implements a semiparametric framework for estimating the intensity function of a spatial point process, where the effect of a target covariate is estimated efficiently in the presence of nonparametric nuisance functions. The method uses a cross-fitting / partialling-out approach to obtain root-$n$ consistent and asymptotically normal estimates of the target parameter.

The code reproduces three sets of numerical results in the paper:

1. **Simulation studies** — Poisson processes, log-Gaussian Cox processes (LGCP), and model misspecification comparisons
2. **Forest ecology analysis** — Beilschmiedia and Capparis tree species in the BCI forest plot
3. **PFAS contamination analysis** — national-level intensity estimation and prediction for per- and polyfluoroalkyl substance (PFAS) measurements across the United States

---

## Repository Structure

```
.
├── code/
│   ├── simulation/
│   │   ├── semi_spp_functions.R              # Core estimation functions
│   │   ├── simulation_functions.R            # Simulation helpers and covariate generators
│   │   ├── generate_gaussian_process.R       # Generate simulated Gaussian covariate fields
│   │   ├── Poisson_simulation.R              # Table 2: Poisson process simulation study
│   │   ├── LGCP_simulation.R                 # Table 3: LGCP simulation study
│   │   └── model_misspecification_simulation.R  # Table 4: model misspecification study
│   │
│   ├── forest_analysis/
│   │   ├── forest_analysis.R                 # Table 5: inference for BCI forest species
│   │   ├── Figure_1.R                        # Figure 1: intensity estimate maps
│   │   └── Figure_2.R                        # Figure 2: estimated nuisance function plot
│   │
│   └── PFAS_analysis/
│       ├── 01_data_prep.R                    # Load and filter PFAS data, build spatial objects
│       ├── 02_covariate_images.R             # Build IDW pixel images for covariates (run once)
│       ├── 03_fit_ipp.R                      # Fit parametric and semiparametric IPP models
│       ├── 04_inference.R                    # SE estimation, inference table, nuisance plot
│       ├── 05_predict_and_plot.R             # National intensity prediction and quantile maps
│       └── 06_zoom_maps.R                    # State-level zoom maps (Alabama, California)
│
├── data/
│   ├── bci.tree1.rdata                       # BCI forest tree census data
│   ├── simulated_Gaussian_field/             # Pre-generated Gaussian covariate fields
│   └── PFAS/
│       ├── IPP Ready Model Data.csv          # PFAS measurements for model fitting
│       └── IPP Ready Extrap Data.csv         # PFAS covariate grid for extrapolation
│
└── output/
    ├── simulation/                           # Simulation result tables (.csv)
    ├── forest/                               # Forest analysis results and figures (.csv, .png)
    └── PFAS/                                 # PFAS fitted objects, maps, and inference (.RData, .png, .csv)
```

---

## Reproducing Results

All scripts assume the working directory is the repository root.

### Simulation Studies (Tables 2–4)

The covariate fields used in all simulations are pre-generated. To regenerate them (requires ~16 GB RAM and ~30 minutes):

```r
source("code/simulation/generate_gaussian_process.R")
```

Then run each simulation study independently:

```r
source("code/simulation/Poisson_simulation.R")              # Table 2 (~3.4 hours)
source("code/simulation/LGCP_simulation.R")                 # Table 3 (~9.9 hours)
source("code/simulation/model_misspecification_simulation.R") # Table 4 (~3.1 hours)
```

Results are saved to `output/simulation/`.

### Forest Ecology Analysis

```r
source("code/forest_analysis/forest_analysis.R")   # Table 5
source("code/forest_analysis/Figure_1.R")          # Figure 1
source("code/forest_analysis/Figure_2.R")          # Figure 2
```

Results are saved to `output/forest/`.

### PFAS Analysis

Run the PFAS pipeline in order. Step 2 is slow (IDW interpolation) and only needs to be run once; its output is reloaded by subsequent steps.

```r
source("code/PFAS_analysis/01_data_prep.R")         # prepare spatial objects
source("code/PFAS_analysis/02_covariate_images.R")  # build covariate images (run once)
source("code/PFAS_analysis/04_inference.R")         # inference (sources 03 automatically)
source("code/PFAS_analysis/05_predict_and_plot.R")  # national maps (sources 03 automatically)
source("code/PFAS_analysis/06_zoom_maps.R")         # state-level zoom maps
```

Results are saved to `output/PFAS/`.

---

## Dependencies

All code is written in R. The following packages are required:

| Package | Use |
|---|---|
| `spatstat` | Point process modelling (`ppm`, `kppm`) |
| `mgcv` | GAM fitting within semiparametric models |
| `gratia` | GAM visualisation |
| `sf` | Spatial data handling |
| `gstat` | IDW spatial interpolation |
| `rnaturalearth`, `rnaturalearthhires` | US state boundary data |
| `ggplot2`, `tidyverse` | Plotting and data manipulation |
| `spam` | Sparse matrix operations for Gaussian field generation |

Install all at once:

```r
install.packages(c("spatstat", "mgcv", "gratia", "sf", "gstat",
                   "rnaturalearth", "rnaturalearthdata", "ggplot2",
                   "tidyverse", "spam"))
remotes::install_github("ropensci/rnaturalearthhires")
```

---

## Citation

If you use this code, please cite:

```
@article{lin2024efficient,
  title={Efficient, Cross-Fitting Estimation of Semiparametric Spatial Point Processes},
  author={Lin, Xindi and others},
  journal={arXiv preprint arXiv:2410.04359},
  year={2024}
}
```
