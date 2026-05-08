# nestprepper <img src="man/figures/logo.png" align="right" height="139" />

[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
The goal of `nestprepper` is to streamline the transition from messy, real-world sea turtle nesting counts to robust biological abundance estimates and Bayesian trend analyses. This package is developed to support NOAA Open Science initiatives by providing a transparent, reproducible pipeline for population status assessments.

## NOAA Disclaimer
This repository is a software product and is not official communication of the National Oceanic and Atmospheric Administration (NOAA), or the United States Department of Commerce. All NOAA GitHub project code is provided "as is" and the user assumes the entire risk as to its quality and performance.

## Overview
`nestprepper` provides a structured workflow to:
1.  **QAQC**: Interactively flag outliers and handle missing monitoring years.
2.  **Biological Conversion**: Transform nest counts into Annual Nesters and Total Adult Females.
3.  **Bayesian Modeling**: Fit State-Space models to estimate population trends ($U$) and process variance ($Q$).
4.  **Visualization**: Generate "Ghost Plots" that communicate population status and uncertainty.

## Data Management & Privacy
To comply with sensitive species data protections, this package does not contain raw GPS or nesting data for specific DPS (Distinct Population Segments). Users should load their own data locally. Simulated data for testing is provided via `simulate_turtle_data()`.

## Installation
You can install the development version of `nestprepper` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("AnnaOrtega-NOAA/nestprepper")
```

## Quick Start
This is a basic example showing how to run a full analysis pipeline:

```r
library(nestprepper)

# 1. Prepare/Clean your data
# This triggers interactive QAQC prompts
clean_data <- prep_nesting_data(my_raw_data)

# 2. Convert Nests to Females
# Based on clutch frequency and remigration intervals
abundance <- calculate_abundance(clean_data, clutch_freq = 5.5, remig_int = 3.0)

# 3. Model Trends (Bayesian State-Space)
results <- run_turtle_model(abundance, iter = 10000, parallel = FALSE)

# 4. Visualize Status
plot_turtle_status(results, abundance, species_name = "Leatherback")
```

## Documentation
For a detailed walkthrough on using your own data and interpreting the mathematical outputs, please see the online tutorial:
`vignette("user-guide", package = "nestprepper")`

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation
If you use `nestprepper` for a DPS assessment or publication, please cite:
Ortega, A. (2026). nestprepper: An R package for sea turtle nesting data integration and trend modeling.
