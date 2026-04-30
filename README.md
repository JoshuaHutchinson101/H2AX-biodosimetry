# γ-H2AX Biodosimetry: Meta-Analysis and Triage Tool

BSc Mathematics and Statistics dissertation, Durham University 2025/26.

This repository contains the key R code that produces tables, figures,
and analyses in the dissertation *"Meta-Analysis of γ-H2AX Decay:
From Hierarchical Pooling to Continuous Temporal Modelling"*, supervised
by Professor Jochen Einbeck.

## Project summary

A complete statistical pipeline for γ-H2AX biodosimetry:

1. **Meta-analytic pooling** of calibration parameters from six
   published studies (46 laboratory-level observations) using a
   multivariate random-effects model with laboratory-level clustering.
2. **Continuous-time decay modelling** via a bi-exponential function
   for β(t) and a single-exponential with floor for α(t), selected
   under formal model comparison against a GAM alternative.
3. **Three-source Monte Carlo dose inversion** combining calibration,
   measurement, and timing uncertainty by refitting the decay models
   to perturbed pooled estimates at each draw.
4. An interactive **Shiny application** that converts a single
   time-stamped foci count into a probabilistic triage classification.

## Repository structure

Meta-analysis and pooling:

- `FE_meta_analysis.R` — Fixed-effect inverse-variance pooling
- `RE_meta_analysis.R` — DerSimonian–Laird random-effects pooling
- `MV_meta_analysis.R` — Multivariate random-effects with lab clustering
- `Cochrans_Q_test.R` — Heterogeneity testing (Cochran's Q, I²)

Continuous-time decay modelling:

- `GAM_fit.R` — Generalised additive model fitting and diagnostics
- `biexp_fit.R` — Bi-exponential model selection (AIC/BIC/LRT)

Dose estimation pipeline:

- `dose_1_setup.R` — Joint covariance matrix construction and Monte Carlo function
- `dose_2_diagnostics.R` — Point estimate confirmation, g-statistic check, model-only distribution
- `dose_3_decomposition.R` — Three-source uncertainty decomposition
- `dose_4_triage.R` — Triage scenario table

Application:

- `Triage_adviser_app.R` — Interactive Shiny application

Scripts are designed to chain — each later script automatically
sources its dependencies if they aren't already in the workspace.

## Reproducing the analysis

Required packages:

    install.packages(c(
      "metafor", "mgcv", "MASS", "Matrix", "shiny",
      "ggplot2", "dplyr", "tidyr", "tibble", "purrr",
      "readxl", "janitor"
    ))

Then run scripts in order.
Tested on R 4.3+ on Windows.

## Live application

The Shiny triage adviser is hosted at:

[link to be added once deployed]

The application takes a measured foci count and elapsed time since
exposure, and returns a point dose estimate, 95% highest-density
interval, and probability vector over Low / Moderate / High triage
categories. Cell count and timing uncertainty are user-adjustable
to reflect operational conditions.

## Author

Joshua Hutchinson  
BSc Mathematics and Statistics, Durham University  
Industrial supervisor: Hannah Mancey, UK Health Security Agency  
Academic supervisor: Professor Jochen Einbeck
